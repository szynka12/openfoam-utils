dimensionSet dimFilter(0, 1, 0, 0, 0, 0, 0);

Info<< "Reading field filter\n" << endl;
volScalarField filter
(
    IOobject
    (
        "filter",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(dimFilter, Zero)
);

if (!filter.headerOk())
{
    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );
    filter = dimensionedScalar("filter", dimFilter, transportProperties);
}

IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

dimensionedScalar viscosity("nu", dimViscosity, transportProperties);

volScalarField pPrim 
(
    IOobject
    (
        "pPrim",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(dimPressure/dimDensity, Zero),
    p.boundaryField().types()
);

volVectorField UPrim
(
    IOobject
    (
        "UPrim",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volTensorField UU
(
    IOobject
    (
        "UU",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
UU = U * U;

volScalarField kStar
(
    IOobject
    (
        "kStar",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar(dimVelocity*dimVelocity, Zero)
 );
