/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2016 OpenFOAM Foundation
    Copyright (C) 2016-2020 OpenCFD Ltd.
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

Class
    Foam::functionObjects::M

Group
    grpFieldFunctionObjects

SourceFiles
    M.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_M_H
#define functionObjects_M_H

#include "fieldExpression.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                          Class M Declaration
\*---------------------------------------------------------------------------*/

class M
:
    public fieldExpression
{
    // Private Member Functions

        //- Calculate the M field and return true if successful
        virtual bool calc();
        const word   k_name;
        const word   k_resolved_name;


public:

    //- Runtime type information
    TypeName("M");


    // Constructors

        //- Construct from Time and dictionary
        M
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- No copy construct
        M(const M&) = delete;

        //- No copy assignment
        void operator=(const M&) = delete;


    //- Destructor
    virtual ~M() = default;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
