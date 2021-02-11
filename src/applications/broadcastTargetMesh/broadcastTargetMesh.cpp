#include "fileOperation.H"
#include "argList.H"
#include "IOobject.H"
#include "IOdictionary.H"
#include "Time.H"



int main(int argc, char *argv[])
{
    
    Foam::argList::addNote
    (
        "Transient solver for incompressible, turbulent flow"
        " of Newtonian fluids on a moving mesh."
    );
    Foam::argList::noParallel();
  
    #include "setRootCase.H"
    #include "createTime.H"
    
    Foam::IOdictionary decomposeParDict 
    (
        Foam::IOobject
        (
            "decomposeParDict",
            runTime.system(),
            runTime,
            Foam::IOobject::MUST_READ_IF_MODIFIED,
            Foam::IOobject::NO_WRITE
        )
    );
    
    Foam::word target_mesh_name 
      = decomposeParDict.getOrDefault<Foam::word>("targetMesh", "targetMesh");

    Foam::fileName target_mesh_path = runTime.rootPath() / runTime.caseName()
        / runTime.constant() / target_mesh_name ;

    Foam::label n_proc 
      = decomposeParDict.get<Foam::label>("numberOfSubdomains");

    for(Foam::label proci = 0; proci < n_proc; proci++)
    {
      Foam::fileName copied_mesh_path = runTime.rootPath() / runTime.caseName()
        / ("processor" + std::to_string(proci)) 
        / runTime.constant() / target_mesh_name;

      Foam::Info << "Copying mesh from: " << Foam::nl 
                 << target_mesh_path << Foam::nl
                 << "to: " << Foam::nl 
                 << copied_mesh_path << Foam::nl << Foam::endl;

      Foam::fileHandler().cp(target_mesh_path, copied_mesh_path);
        
    }


    return 0;
};
