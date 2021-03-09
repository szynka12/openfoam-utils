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
    Foam::functionObjects::visRatio

Group
    grpFieldFunctionObjects

SourceFiles
    visRatio.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_visRatio_H
#define functionObjects_visRatio_H

#include "fieldExpression.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                          Class visRatio Declaration
\*---------------------------------------------------------------------------*/

class visRatio
:
    public fieldExpression
{
    // Private Member Functions

        //- Calculate the visRatio field and return true if successful
        virtual bool calc();


public:

    //- Runtime type information
    TypeName("visRatio");


    // Constructors

        //- Construct from Time and dictionary
        visRatio
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- No copy construct
        visRatio(const visRatio&) = delete;

        //- No copy assignment
        void operator=(const visRatio&) = delete;


    //- Destructor
    virtual ~visRatio() = default;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
