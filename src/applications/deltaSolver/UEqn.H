// Solve the Momentum equation

MRF.correctBoundaryVelocity(U);

kStar = 0.5 * UPrim & UPrim;

tmp<volTensorField> tgradU = fvc::grad(U);
const volTensorField& gradU = tgradU();
 
tmp<volScalarField> Smag
(
    dev(twoSymm(gradU)) && gradU
);

tmp<volSymmTensorField> S(symm(gradU));



tmp<fvVectorMatrix> tUEqn
(
    fvm::ddt(U) + fvm::div(phi, U)
  + MRF.DDt(U)
  - fvm::laplacian(viscosity, U)
 ==
    fvOptions(U)
  - fvc::div(UU - U*U)  //leoneard term
  - fvc::div(   
                UPrim*U + U*UPrim + UPrim*UPrim 
              - S * 2 * dimensionedScalar(dimless/dimTime, 0.1)
                / (Smag + dimensionedScalar(Smag()().dimensions(), SMALL)) 
                * kStar  
            )
);
fvVectorMatrix& UEqn = tUEqn.ref();

UEqn.relax();

fvOptions.constrain(UEqn);

if (pimple.momentumPredictor())
{
    solve(UEqn == -fvc::grad(p));

    fvOptions.correct(U);
}
