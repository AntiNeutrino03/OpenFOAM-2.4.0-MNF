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

\*---------------------------------------------------------------------------*/

#include "radialStripToCell.H"
#include "polyMesh.H"

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(radialStripToCell, 0);

addToRunTimeSelectionTable(topoSetSource, radialStripToCell, word);

addToRunTimeSelectionTable(topoSetSource, radialStripToCell, istream);

}


Foam::topoSetSource::addToUsageTable Foam::radialStripToCell::usage_
(
    radialStripToCell::typeName,
    "\n    Usage: radialStripToCell (centreX centreY centreZ) \n\n"
    "    outerRadius innerRadius \n\n"
    "    Select all cells with cellCentre within bounding region\n\n"
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::radialStripToCell::combine(topoSet& set, const bool add) const
{
    const pointField& ctrs = mesh_.cellCentres();

    forAll(ctrs, cellI)
    {
        scalar offset = sqrt(
                                sqr(centre_.y() - ctrs[cellI].y()) +
                                sqr(centre_.z() - ctrs[cellI].z())
                            );
        
        scalar cellAngle = atan(
                                    (centre_.z() - ctrs[cellI].z()) /
                                    (centre_.y() - ctrs[cellI].y())
                                )*57.2958;
                                
        vector offsets = ctrs[cellI] - centre_;
        
//         if(offsets.z() <= 0 && offsets.y() >= 0)
//         {
//             Info << "offsets = " << offsets << endl;
//             Info << "cellAngle = " << cellAngle << endl;
//         }
                                
        if(offsets.z() >= 0)
        {
            if(offsets.y() >= 0) //1st quadrant
            {
                //ok
            }
            else //2nd quadrant
            {
                cellAngle = 180 + cellAngle; //ok
            }
        }
        else
        {
           if(offsets.y() <= 0) //3rd quadrant
            {
                cellAngle = 180 + cellAngle; //ok
            }
            else //4th quadrant
            {
                cellAngle = 360 + cellAngle; //ok
            } 
        }
        
        if(offset <= outerRadius_ && offset >= innerRadius_)
        {
            if(ctrs[cellI].x() >= axialStart_ && ctrs[cellI].x() <= axialEnd_)
            {
                if(cellAngle >= angleStart_ && cellAngle <= angleEnd_)
                {
                    addOrDelete(set, cellI, add);
                }
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::radialStripToCell::radialStripToCell
(
    const polyMesh& mesh,
    const vector& centre,
    const scalar outerRadius,
    const scalar innerRadius,
    const scalar axialStart,
    const scalar axialEnd,
    const scalar angleStart,
    const scalar angleEnd
)
:
    topoSetSource(mesh),
    centre_(centre),
    outerRadius_(outerRadius),
    innerRadius_(innerRadius),
    axialStart_(axialStart),
    axialEnd_(axialEnd),
    angleStart_(angleStart),
    angleEnd_(angleEnd)
{}


// Construct from dictionary
Foam::radialStripToCell::radialStripToCell
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    topoSetSource(mesh),
    centre_(dict.lookup("centre")),
    outerRadius_(readScalar(dict.lookup("outerRadius"))),
    innerRadius_(readScalar(dict.lookup("innerRadius"))),
    axialStart_(readScalar(dict.lookup("axialStart"))),
    axialEnd_(readScalar(dict.lookup("axialEnd"))),
    angleStart_(readScalar(dict.lookup("angularStart"))),
    angleEnd_(readScalar(dict.lookup("angularEnd")))
{}


// Construct from Istream
Foam::radialStripToCell::radialStripToCell
(
    const polyMesh& mesh,
    Istream& is
)
:
    topoSetSource(mesh),
    centre_(checkIs(is)),
    outerRadius_(readScalar(checkIs(is))),
    innerRadius_(readScalar(checkIs(is))),
    axialStart_(readScalar(checkIs(is))),
    axialEnd_(readScalar(checkIs(is))),
    angleStart_(readScalar(checkIs(is))),
    angleEnd_(readScalar(checkIs(is)))
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radialStripToCell::~radialStripToCell()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::radialStripToCell::applyToSet
(
    const topoSetSource::setAction action,
    topoSet& set
) const
{
    if ((action == topoSetSource::NEW) || (action == topoSetSource::ADD))
    {
        Info<< "    Adding cells within defined region " << endl;

        combine(set, true);
    }
    else if (action == topoSetSource::DELETE)
    {
        Info<< "    Removing cells within defined region " << endl;

        combine(set, false);
    }
}


// ************************************************************************* //
