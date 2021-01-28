#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "CorrectPhi.H"
#include "fvOptions.H"

template <typename type>
void laplace_filter_(GeometricField< type, fvPatchField, volMesh >& field,
    const volScalarField& filter)
{
  
  IOdictionary fvSolution
  (
      IOobject
      (
          "fvSolution",
          field.mesh().time().system(),
          field().mesh(),
          IOobject::MUST_READ_IF_MODIFIED,
          IOobject::NO_WRITE
      )
  );
  

  fvMatrix<type> eq
  (
     - fvm::laplacian(filter, field)
     ==
     - fvm::Sp(dimensionedScalar(dimless, 1.0), field)
     + field 
  );
  eq.solve(fvSolution.subDict("solvers").subDict("filter"));
}

template <typename type>
void sph_filter(GeometricField< type, fvPatchField, volMesh >& field,
    const volScalarField& filter)
{
  laplace_filter_(field, sqr(filter)/40);
}

//tmp<volTensorField> leonard_stress(const )

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Transient solver for incompressible, turbulent flow"
        " of Newtonian fluids on a moving mesh."
    );
    // clang-format off
    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "initContinuityErrs.H"
    #include "createDyMControls.H"

    #include "createFields.H"
    
    // explict LES fields
    #include "createLESFields.H"
    
    #include "createUfIfPresent.H"
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"
    // clang-format on

    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        // clang-format off
        #include "readDyMControls.H"
        #include "CourantNo.H"
        #include "setDeltaT.H"
        // clang-format on
        
        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            if (pimple.firstIter() || moveMeshOuterCorrectors)
            {
                // Do any mesh changes
                mesh.controlledUpdate();

                if (mesh.changing())
                {
                    MRF.update();

                    if (correctPhi)
                    {
                        // Calculate absolute flux
                        // from the mapped surface velocity
                        phi = mesh.Sf() & Uf();

                        #include "correctPhi.H"

                        // Make the flux relative to the mesh motion
                        fvc::makeRelative(phi, U);
                    }

                    if (checkMeshCourantNo)
                    {
                        #include "meshCourantNo.H"
                    }
                }
            }

            #include "UEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                laminarTransport.correct();
                turbulence->correct();
            }

            
        }

        pPrim = p;
        sph_filter<scalar>(p, filter);
        pPrim -= p;

        UPrim = U;
        sph_filter<vector>(U, filter);
        UPrim -= U;
        sph_filter<vector>(UPrim, filter);

        //update UU for next iteration
        UU = U * U;
        sph_filter<tensor>(UU, filter);


        phi = fvc::flux(U);
        #include "continuityErrs.H"

        runTime.write();

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
