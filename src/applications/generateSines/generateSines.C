#include "fvCFD.H"
#include "fvOptions.H"
#include "simpleControl.H"


struct Waves
{
  const List<scalar> f;
  const List<scalar> As;
  const List<scalar> Ac;

  Waves(IOdictionary& dict)
    :
  f(dict.get<List<scalar>>("f")),
  As(dict.get<List<scalar>>("sin")),
  Ac(dict.get<List<scalar>>("cos"))
  {
    if (f.size() != As.size() || f.size() != Ac.size())
    {
      FatalError << "All lists in generateSinesDict must have the same size!"
                 << endl;
    }
  }
  
  scalar eval(const vector& x)
  { 
    scalar value = 0;

    const scalar p = Foam::mag(x);

    forAll(f, i)
    {
      value 
        += As[i] * Foam::sin(2 * Foam::constant::mathematical::pi * f[i] * p)
        + Ac[i] * Foam::cos(2 * Foam::constant::mathematical::pi * f[i]* p);
    }
    return value;
  }

};




// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Generate sine waves according to generateSinesDict"
    );
 
    #include "postProcess.H"
 
    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
 
    simpleControl simple(mesh);

    #include "createFields.H"

    IOdictionary sinesDict
    (
        IOobject
        (
            "generateSinesDict",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );
    
    Waves waves(sinesDict);
 
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
    Info<< "\nGenerating sine waves...\n" << endl;
 
    while (simple.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;
        
        forAll(p, i)
        {
          p[i] = waves.eval(mesh.cellCentres()[i]);
        }
        


        if (runTime.writeTime()) { runTime.write(); }

        runTime.printExecutionTime(Info);
    }
 
    Info<< "End\n" << endl;
 
    return 0;
}
 
 
// ************************************************************************* //
