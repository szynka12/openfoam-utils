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


    // For each time step read all fields
    forAll(timeDirs, timeI)
    {

        runTime.setTime(timeDirs[timeI], timeI);
        
        Info << "Reading fields for time " << runTime.timeName() << endl;

        // if (timeI==0){
        //     Info << "Skipping " << runTime.timeName() << endl;
        //     continue;
        //     }


        #include "readFields.H"
        
        #include "calculateFields.H"

        // Average fields over channel down to a line
        Info<< "Collapsing fields for time " << runTime.timeName() << endl;
        #include "collapse.H"
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
