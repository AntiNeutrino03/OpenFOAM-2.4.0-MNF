/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
     \\/     M anipulation  |
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
    Foam::radialStripToCell

Description
    A topoSetSource to select cells based on cell centres inside a radial range
    and a circumferencial range.
    
Usage
    \verbatim
    {
        name    strip1;
        type    cellSet;
        action  new;
        source  radialStripToCell;
        sourceInfo
        {
            centre (0.1802 -0.01640 0);
            outerRadius 0.1200;
            innerRadius 0.10555;
            axialStart 0.1800;
            axialEnd 0.1805;
            angularStart 160;
            angularEnd 200;
        }
    }
    \endverbatim

SourceFiles
    radialStripToCell.C

\*---------------------------------------------------------------------------*/

#ifndef radialStripToCell_H
#define radialStripToCell_H

#include "topoSetSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class radialStripToCell Declaration
\*---------------------------------------------------------------------------*/

class radialStripToCell
:
    public topoSetSource
{

    // Private data

        //- Add usage string
        static addToUsageTable usage_;


        //- Centre
        vector centre_;
        
        //- Outer radius
        scalar outerRadius_;
        
        //- Inner radius
        scalar innerRadius_;
        
        //- Axial start
        scalar axialStart_;
        
        //- Axial end
        scalar axialEnd_;
        
        //- Angular start
        scalar angleStart_;
        
        //- Angular end
        scalar angleEnd_;


    // Private Member Functions

        void combine(topoSet& set, const bool add) const;


public:

    //- Runtime type information
    TypeName("radialStripToCell");

    // Constructors

        //- Construct from components
        radialStripToCell
        (
            const polyMesh& mesh,
            const vector& centre,
            const scalar outerRadius,
            const scalar innerRadius,
            const scalar axialStart,
            const scalar axialEnd,
            const scalar angleStart,
            const scalar angleEnd
        );

        //- Construct from dictionary
        radialStripToCell
        (
            const polyMesh& mesh,
            const dictionary& dict
        );

        //- Construct from Istream
        radialStripToCell
        (
            const polyMesh& mesh,
            Istream&
        );


    //- Destructor
    virtual ~radialStripToCell();


    // Member Functions

        virtual sourceType setType() const
        {
            return CELLSETSOURCE;
        }

        virtual void applyToSet
        (
            const topoSetSource::setAction action,
            topoSet&
        ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
