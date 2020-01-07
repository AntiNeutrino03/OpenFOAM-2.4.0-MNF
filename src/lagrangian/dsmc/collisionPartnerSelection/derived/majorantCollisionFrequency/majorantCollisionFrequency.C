/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2005 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    majorantCollisionFrequency

Description

\*----------------------------------------------------------------------------*/

#include "majorantCollisionFrequency.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(majorantCollisionFrequency, 0);

addToRunTimeSelectionTable
(collisionPartnerSelection, majorantCollisionFrequency, dictionary);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from components
majorantCollisionFrequency::majorantCollisionFrequency
(
    const polyMesh& mesh,
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    collisionPartnerSelection(mesh, cloud, dict),
    infoCounter_(0)
//     propsDict_(dict.subDict(typeName + "Properties"))
{}



// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

majorantCollisionFrequency::~majorantCollisionFrequency()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void majorantCollisionFrequency::initialConfiguration()
{

}

void majorantCollisionFrequency::collide()
{
    if (!cloud_.binaryCollision().active())
    {
        return;
    }

    // Temporary storage for subCells
    List<DynamicList<label> > subCells(8);

    const scalar& deltaT = cloud_.mesh().time().deltaTValue();
    label collisions = 0;

    const List<DynamicList<dsmcParcel*> >& cellOccupancy = 
                                                cloud_.cellOccupancy();

    const polyMesh& mesh = cloud_.mesh();

    forAll(cellOccupancy, cellI)
    {
        const DynamicList<dsmcParcel*>& cellParcels(cellOccupancy[cellI]);
        const scalar& cellVolume = mesh.cellVolumes()[cellI];
        scalar sumLocalTimeStep = 0;
        scalar nuMax = 0.0;
        scalar localTimeStep = 0.0;

        label nC(cellParcels.size());

        if (nC > 1)
        {
            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            // Assign particles to one of 8 Cartesian subCells

            // Clear temporary lists
            forAll(subCells, i)
            {
                subCells[i].clear();
            }

            // Inverse addressing specifying which subCell a parcel is in
            List<label> whichSubCell(cellParcels.size());

            point cC = mesh.cellCentres()[cellI];

            forAll(cellParcels, i)
            {
                const dsmcParcel& p = *cellParcels[i];

                vector relPos = p.position() - cC;

                label subCell =
                    pos(relPos.x()) + 2*pos(relPos.y()) + 4*pos(relPos.z());

                subCells[subCell].append(i);

                whichSubCell[i] = subCell;
            }

            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            while(1 < 2)
            {

                scalar sigmaTcRMax = cloud_.sigmaTcRMax()[cellI];

                scalar R = -1.0;

                if(cloud_.axisymmetric())
                {
                    scalar RWF = 0.0;

                    const point& cC = mesh.cellCentres()[cellI];

                    scalar radius = cC.y();

                    RWF = 1.0 + cloud_.maxRWF()
                                            *(radius/cloud_.radialExtent());

                    nuMax = 0.5*nC*(nC - 1)*cloud_.nParticle()
                                *RWF*sigmaTcRMax/cellVolume;
                }
                else
                {
                    nuMax = 0.5*nC*(nC - 1)*cloud_.nParticle()
                                    *sigmaTcRMax/cellVolume;
                }
                
                do
                {
                    R = rndGen_.scalar01();
                } while(R < VSMALL);

                if(nuMax > VSMALL)
                {
                    localTimeStep = -log(R)/nuMax;
                }
                
                sumLocalTimeStep += localTimeStep;
                
                if(sumLocalTimeStep >= deltaT)
                {
                    break;
                }
            
                // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                // subCell candidate selection procedure

                // Select the first collision candidate
                label candidateP = rndGen_.integer(0, nC - 1);

                // Declare the second collision candidate
                label candidateQ = -1;

                const List<label>& subCellPs = 
                                        subCells[whichSubCell[candidateP]];

                label nSC = subCellPs.size();

                if (nSC > 1)
				{
                    // If there are two or more particle in a subCell, choose
                    // another from the same cell.  If the same candidate is
                    // chosen, choose again.

                    do
                    {
                        candidateQ = subCellPs[rndGen_.integer(0, nSC - 1)];

                    } while (candidateP == candidateQ);
                }
                else
                {
                    // Select a possible second collision candidate from the
                    // whole cell.  If the same candidate is chosen, choose
                    // again.

                    do
                    {
                        candidateQ = rndGen_.integer(0, nC - 1);

                    } while (candidateP == candidateQ);
                }

                dsmcParcel& parcelP = *cellParcels[candidateP];
                dsmcParcel& parcelQ = *cellParcels[candidateQ];

                label chargeP = -2;
                label chargeQ = -2;

                chargeP = cloud_.constProps(parcelP.typeId()).charge();
                chargeQ = cloud_.constProps(parcelQ.typeId()).charge();
                
                //do not allow electron-electron collisions
                
                if(!(chargeP == -1 && chargeQ == -1))
                {    
                    scalar sigmaTcR = cloud_.binaryCollision().sigmaTcR
                    (
                        parcelP,
                        parcelQ
                    );


                    if (sigmaTcR > cloud_.sigmaTcRMax()[cellI])
                    {
                        cloud_.sigmaTcRMax()[cellI] = sigmaTcR;
                    }

                    if (sigmaTcR/sigmaTcRMax > rndGen_.scalar01())
                    {
                        // chemical reactions

                        // find which reaction model parcel p and q should use
                        label rMId = cloud_.reactions().returnModelId
                                                        (parcelP, parcelQ);

                        if(rMId != -1)
                        {
                                cloud_.reactions().reactions()[rMId]->reaction
                                (
                                    parcelP,
                                    parcelQ
                                );                                   
                            // if reaction unsuccessful use conventional
                            //    collision model
                            if(cloud_.reactions().reactions()[rMId]->relax())
                            {
                                cloud_.binaryCollision().collide
                                (
                                    parcelP,
                                    parcelQ,
                                    cellI
                                );
                            }
                           
                        }
                        
                        // if reaction model not found, use conventional 
                        // collision model
                        else 
                        {
                            cloud_.binaryCollision().collide
                            (
                                parcelP,
                                parcelQ,
                                cellI
                            );
                        }

                        collisions++;
                    }
                    //otherwise a fictitious collision
                }
            } 
        }
    }

    reduce(collisions, sumOp<label>());
    cloud_.sigmaTcRMax().correctBoundaryConditions();
    
    infoCounter_++;
        
    //if (sumlocalTimeStep <= deltaT)
   if(infoCounter_ >= cloud_.nTerminalOutputs())
   {
            Info<< "    Collisions                      = "
                << collisions << nl
                << endl;
                
            infoCounter_ = 0; 
   }
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
