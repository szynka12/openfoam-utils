/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    postChannel

Group
    grpPostProcessingUtilities

Description
    Post-process data from channel flow calculations.


\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "channelIndex.H"
#include "makeGraph.H"

#include "OSspecific.H"

template <class VolField>
inline bool check_header(IOobject& head)
{
    bool result = !head.typeHeaderOk<VolField>();
    if (result){ Info << "    " << "No " << head.name() << " field" << endl; }
    return result;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Post-process data from channel flow calculations"
    );

    argList::noParallel();
    timeSelector::addOptions();

    #include "setRootCase.H"
    #include "createTime.H"

    // Get times list
    instantList timeDirs = timeSelector::select0(runTime, args);

    #include "createNamedMesh.H"
    #include "readTransportProperties.H"

    dimensionedScalar nu("nu", dimViscosity, transportProperties);

    const word& gFormat = runTime.graphFormat();

    // Setup channel indexing for averaging over channel down to a line

    IOdictionary channelDict
    (
        IOobject
        (
            "postChannelDict",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );
    channelIndex channelIndexing(mesh, channelDict);

    const word U_name = channelDict.getOrDefault<word>("U", "UMean");
    const word production_name = channelDict.getOrDefault<word>(
        "production", "production");

    const word dissipation_name = channelDict.getOrDefault<word>(
        "dissipation", "dissipation");
    
    const word vpg_name = channelDict.getOrDefault<word>(
        "velocityPressureGradient", "velocityPressureGradient");
    
    const word mD_name = channelDict.getOrDefault<word>(
        "molecularDiffusion", "molecularDiffusion");

    const word tD_name = channelDict.getOrDefault<word>(
        "turbulentDiffusion", "turbulentDiffusion");

    bool only_velocity = channelDict.getOrDefault<bool>("onlyVel", false);
    
    // For each time step read all fields
    forAll(timeDirs, timeI)
    {

        runTime.setTime(timeDirs[timeI], timeI);
        
        Info << "Reading fields for time " << runTime.timeName() << endl;

        IOobject UMeanHeader(U_name, runTime.timeName(), mesh, IOobject::MUST_READ);
        if (check_header<volVectorField>(UMeanHeader)) {continue;}
        volVectorField UMean(UMeanHeader, mesh);

        IOobject productionHeader(
            production_name, runTime.timeName(), mesh, IOobject::MUST_READ);
        if (check_header<volScalarField>(productionHeader)) {continue;}
        volScalarField production(productionHeader, mesh);

        IOobject dissipationHeader(
            dissipation_name, runTime.timeName(), mesh, IOobject::MUST_READ);
        if (check_header<volScalarField>(dissipationHeader)) {continue;}
        volScalarField dissipation(dissipationHeader, mesh);

        IOobject vpgHeader(
            vpg_name, runTime.timeName(), mesh, IOobject::MUST_READ);
        if (check_header<volScalarField>(vpgHeader)) {continue;}
        volScalarField vpg(vpgHeader, mesh);

        IOobject molDiffHeader(
            mD_name, runTime.timeName(), mesh, IOobject::MUST_READ);
        if (check_header<volScalarField>(molDiffHeader)) {continue;}
        volScalarField molDiff(molDiffHeader, mesh);

        IOobject turbDiffHeader(
            tD_name, runTime.timeName(), mesh, IOobject::MUST_READ);
        if (check_header<volScalarField>(turbDiffHeader)) {continue;}
        volScalarField turbDiff(turbDiffHeader, mesh);
    

        //#include "readFields.H"
        
        #include "calculateFields.H"

        // Average fields over channel down to a line
        Info<< "Collapsing fields for time " << runTime.timeName() << endl;
        #include "collapse.H"
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
