#include "fileOperation.H"
#include "argList.H"
#include "IOobject.H"
#include "IOdictionary.H"
#include "Time.H"
#include "PtrList.H"
#include "fvMesh.H"
#include "GeometricField.H"
#include "volFieldsFwd.H"
#include "fvCFD.H"
#include "IOobjectList.H"


template <class volFieldType>
void reconstruct_field(IOobjectList& objects, const Time& target_time,
    const fvMesh& mesh, const PtrList<Time>& databases, 
    const volScalarField& volume)
{
  const wordList names  = objects.lookupClass<volFieldType>().names();
  
  forAll(names, namei)
  {

    volFieldType processor_field
    (
        Foam::IOobject
        (
            names[namei],
            databases[0].timeName(),
            databases[0],
            Foam::IOobject::MUST_READ,
            Foam::IOobject::AUTO_WRITE
        ),
        mesh
    );

    for (label proci = 1; proci < databases.size(); proci++)
    {
      processor_field += volFieldType(
          Foam::IOobject
          (
              names[namei],
              databases[proci].timeName(),
              databases[proci],
              Foam::IOobject::MUST_READ,
              Foam::IOobject::AUTO_WRITE
          ),
          mesh
      );
    }
    
    IOobject reconstructed_header 
      = IOobject(names[namei], target_time.timeName(), target_time,
          IOobject::NO_READ, IOobject::AUTO_WRITE);
    
    volFieldType reconstructed_field(reconstructed_header, processor_field);

    reconstructed_field *= (1.0 / (volume + SMALL ));

    reconstructed_field.write();
  }

}



int main(int argc, char *argv[])
{
    
    Foam::argList::addNote
    (
        "Gather results from cellFilter function object."
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
    
    const Foam::fileName target_mesh_name 
      = decomposeParDict.getOrDefault<Foam::word>("targetMesh", "targetMesh");

    const Foam::fileName target_root_path 
      = runTime.rootPath() / runTime.caseName() / runTime.constant();

    const Foam::fileName target_mesh_path =  target_root_path/target_mesh_name;

    Foam::label n_proc 
      = decomposeParDict.get<Foam::label>("numberOfSubdomains");

    
    // create time db for the target_mesh
    Foam::Time target_run_time(
        runTime.rootPath() / runTime.caseName() / runTime.constant(), 
        target_mesh_name);
    
    // Create the processor databases
    Foam::PtrList<Foam::Time> databases(n_proc);
  
    forAll(databases, proci)
    {
        databases.set
        (
            proci,
            new Foam::Time(
                runTime.rootPath() / runTime.caseName() 
                / ("processor" + Foam::name(proci)) / runTime.constant(),
                target_mesh_name )
        );
    }

    // Read targt mesh
    Foam::autoPtr<Foam::fvMesh> meshPtr(nullptr);
    Foam::word regionName = Foam::fvMesh::defaultRegion;
    
    meshPtr.reset
    (
        new Foam::fvMesh
        (
            Foam::IOobject
            (
                regionName,
                target_run_time.timeName(),
                target_run_time,
                Foam::IOobject::MUST_READ
            )
        )
    );
    Foam::fvMesh& mesh = meshPtr();

    // Create the volume field for the target mesh mesh
    Foam::volScalarField volume
    (
        Foam::IOobject
        (
            "FilterVolume",
            target_run_time.timeName(),
            mesh,
            Foam::IOobject::NO_READ,
            Foam::IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dimless, Zero)
    );
    
    Foam::Info << "Summing volumes..." << endl;
    // First sum the volumes in constant directories
    for(Foam::label proci = 0; proci < n_proc; proci++)
    {
      Foam::Info << "    Reading vol for processor " << Foam::name(proci) << endl;
      
      const Foam::volScalarField local_volume
      (
          Foam::IOobject
          (
              "FilterVolume",
              databases[proci].constant(),
              databases[proci],
              Foam::IOobject::MUST_READ,
              Foam::IOobject::NO_WRITE
          ),
          mesh
      );

      volume += local_volume;
    }

    Foam::Info << "Done!" << endl;

    Foam::Info << "Writing volume field ..." << endl;
    volume.write();

    // Get the times from the master processor
    instantList timeDirs = databases[0].times();

    forAll(timeDirs, timei)
    {


      // Set time for global database
      target_run_time.setTime(timeDirs[timei], timei);
      
      forAll(databases, proci)
      {
        databases[proci].setTime(timeDirs[timei], timei);
      }
      if ( target_run_time.timeName() == "constant")
        continue;

      Foam::Info<< "Time = " << target_run_time.timeName() << endl << endl;
      
      Foam::IOobjectList objects( databases[0], databases[0].timeName() );
      
      reconstruct_field<volScalarField>(objects, target_run_time, mesh, 
          databases, volume);

      reconstruct_field<volVectorField>(objects, target_run_time, mesh, 
          databases, volume);
      
      reconstruct_field<volSymmTensorField>(objects, target_run_time, mesh, 
          databases, volume);
  
      reconstruct_field<volTensorField>(objects, target_run_time, mesh, 
          databases, volume);
    }

    return 0;
};
