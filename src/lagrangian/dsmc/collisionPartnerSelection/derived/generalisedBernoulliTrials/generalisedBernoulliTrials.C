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

\*----------------------------------------------------------------------------*/

#include "generalisedBernoulliTrials.H"
#include "addToRunTimeSelectionTable.H"
#include <algorithm>

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(generalisedBernoulliTrials, 0);

addToRunTimeSelectionTable(collisionPartnerSelection, generalisedBernoulliTrials, dictionary);



// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from components
generalisedBernoulliTrials::generalisedBernoulliTrials
(
    const polyMesh& mesh,
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    collisionPartnerSelection(mesh, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    moveToTop_(propsDict_.lookup("sortParticleList")),
    inverseSel_(),
    definedSel_(),
    inverseSelectionFraction_
    (
        propsDict_.lookupOrDefault<scalar>("inverseSelectionFraction", 0.0)
    ),
    nSel_
    (
        propsDict_.lookupOrDefault<scalar>("nParcelsToSelect", 0.0)
    )
{}



// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

generalisedBernoulliTrials::~generalisedBernoulliTrials()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void generalisedBernoulliTrials::initialConfiguration()
{
    bool inverseSel;
    
    if(inverseSelectionFraction_ < VSMALL)
    {
        inverseSel_ = false;
    }
    else
    {
        inverseSel_ = true;
    }
    
    if(nSel_ < VSMALL)
    {
        definedSel_ = false;
    }
    else
    {
        definedSel_ = true;
    }
    
    if(inverseSel_ && definedSel_)
    {
        FatalErrorIn("generalisedBernoulliTrials::initialConfiguration()")
        << "Two methods for selecting the number of particles to sort "
        << "have been selected, there should only be one. Remove either"
        << " 'inverseSelectionFraction' or 'nParcelsToSort' option"
        << nl
        << exit(FatalError);
    }
    
    if(!inverseSel_ && !definedSel_)
    {
        FatalErrorIn("generalisedBernoulliTrials::initialConfiguration()")
        << "No for selecting the number of particles to sort "
        << "has been selected, there must be be one. Add either"
        << " 'inverseSelectionFraction' or 'nParcelsToSort' option"
        << nl
        << exit(FatalError);
    }
}

void generalisedBernoulliTrials::collide()
{
    if (!cloud_.binaryCollision().active())
    {
        return;
    }
    
    // Temporary storage for subCells
    List<DynamicList<label> > subCells(8);

    const scalar& deltaT = cloud_.mesh().time().deltaTValue();

    label collisions = 0;

    const List<DynamicList<dsmcParcel*> >& cellOccupancy = cloud_.cellOccupancy();

    const polyMesh& mesh = cloud_.mesh();
    
    forAll(cellOccupancy, cellI)
    {
        DynamicList<dsmcParcel*> cellParcels(cellOccupancy[cellI]);

        scalar nC(cellParcels.size());
        
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
        
        scalar nSel = 0.0;
        
        if(inverseSel_)
        {
            //number of parcels to select per subCell
            nSel = nC/inverseSelectionFraction_/8.0;
        }
        
        if(definedSel_)
        {
            //number of parcels to select per subCell
            nSel = nSel_/8.0;
        }
    
        if(nC < VSMALL)
        {
            nSel = 1.0;
        }
        
        scalar nParticle = cloud_.nParticle();
        scalar prob1 = (nParticle*deltaT)/(mesh.cellVolumes()[cellI]/8.0);
        label k = -1;
        scalar kPrime = -1;
        label candidateP = -1;
        label candidateQ = -1;
        
        forAll(subCells, i)
        {
            label nCS = subCells[i].size();
            
            if(nSel < (nCS - 1.0))
            {
                //follow the GBT algorithm
                
                if(moveToTop_)
                {
                    DynamicList<label> selectedParcels(0);
                    
                    label candidate = 0;
                    
                    //randomly choose particles to move to top of each subCell list
                    for(label n = 0; n < nSel - 1; n++)
                    {
                        bool found = false;
                        
                        do
                        {
                            candidate = rndGen_.integer(0, nC - 1);
                            
                            //don't allow the same particle to be chosen twice
                            found = (std::find(selectedParcels.begin(), selectedParcels.end(), candidate) != selectedParcels.end());
                            
                        }while(found);
                        
                        selectedParcels.append(candidate);
                    }
                    
                    //move the particles to the top of the list
                    forAll(selectedParcels, j)
                    {
                        dsmcParcel& storeParcel1 = *cellParcels[selectedParcels[j]];
                        dsmcParcel& storeParcel2 = *cellParcels[subCells[i][j]];
                        
                        cellParcels[subCells[i][j]] = &storeParcel1;
                        cellParcels[selectedParcels[j]] = &storeParcel2;
                    }
                }
                
                //run SBT procedure from 0 to nSel-1
                for(label n = 0; n < nSel - 1; n++)
                {
                    scalar RWF = 1.0;
                    
                    if(cloud_.axisymmetric())
                    {
                        const point& cC = mesh.cellCentres()[cellI];
                        
                        scalar radius = cC.y();
                        
                        RWF = 1.0 + cloud_.maxRWF()*(radius/cloud_.radialExtent());
                        
                        nParticle *= RWF;
                    }

                    candidateP = n;
                    kPrime = ((nCS*(nCS-1.0)))/(nSel*(2.0*nCS-nSel-1.0));
                    k = nCS - 1 - n;
                    candidateQ = n + rndGen_.integer(1, k);
                    
                    dsmcParcel& parcelP = *cellParcels[subCells[i][candidateP]];
                    dsmcParcel& parcelQ = *cellParcels[subCells[i][candidateQ]];

                    scalar sigmaTcR = cloud_.binaryCollision().sigmaTcR
                    (
                        parcelP,
                        parcelQ
                    );
                
                    scalar probability = k*kPrime*prob1*sigmaTcR;

                    if (probability > rndGen_.scalar01())
                    {
                        // chemical reactions

                        // find which reaction model parcel p and q should use
                        label rMId = cloud_.reactions().returnModelId(parcelP, parcelQ);

                        if(rMId != -1)
                        {
                            // try to react molecules
                            if(cloud_.reactions().reactions()[rMId]->reactWithLists())
                            {

                            }
                            else
                            {
                                cloud_.reactions().reactions()[rMId]->reaction
                                (
                                    parcelP,
                                    parcelQ
                                );                                    
                            }
                            // if reaction unsuccessful use conventional collision model
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
                        else // if reaction model not found, use conventional collision model
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
                }
            }
            else
            {
                //follow the SBT aglorithm
                if (nC > 1)
                {   
                    label nCS = subCells[i].size();
                    
                    scalar nParticle = cloud_.nParticle();
                    
                    scalar RWF = 1.0;
                    
                    if(cloud_.axisymmetric())
                    {
                        const point& cC = mesh.cellCentres()[cellI];
                        
                        scalar radius = cC.y();
                        
                        RWF = 1.0 + cloud_.maxRWF()*(radius/cloud_.radialExtent());
                        
                        nParticle *= RWF;
                    }
                    
                    scalar prob1 = (nParticle*deltaT)/(mesh.cellVolumes()[cellI]/8.0);
                    label k = -1;
                    label candidateP = -1;
                    label candidateQ = -1;
                    
                    // loop over sub cells
                    forAll(subCells, i)
                    {     
                        if(nCS > 1)
                        {       
                            for(label p = 0 ; p < nCS - 1 ; p++)
                            {                
                                // Select the first collision candidate
                                candidateP = p;
                            
                                k = nC - 1 - p;
                                candidateQ = p + rndGen_.integer(1, k);
                                
                                dsmcParcel& parcelP = *cellParcels[subCells[i][candidateP]];
                                dsmcParcel& parcelQ = *cellParcels[subCells[i][candidateQ]]; 

                                scalar sigmaTcR = cloud_.binaryCollision().sigmaTcR
                                (
                                    parcelP,
                                    parcelQ
                                );
                            
                                scalar probability = k*prob1*sigmaTcR;

                                if (probability > rndGen_.scalar01())
                                {
                                    // chemical reactions

                                    // find which reaction model parcel p and q should use
                                    label rMId = cloud_.reactions().returnModelId(parcelP, parcelQ);

                                    if(rMId != -1)
                                    {
                                        // try to react molecules
                                        if(cloud_.reactions().reactions()[rMId]->reactWithLists())
                                        {

                                        }
                                        else
                                        {
                                            cloud_.reactions().reactions()[rMId]->reaction
                                            (
                                                parcelP,
                                                parcelQ
                                            );                                    
                                        }
                                        // if reaction unsuccessful use conventional collision model
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
                                    else // if reaction model not found, use conventional collision model
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
                            }
                        }
                    }
                }
            }
        }
    }

    reduce(collisions, sumOp<label>());

//     sigmaTcRMax_.correctBoundaryConditions();

    if (collisions > 0)
    {
        Info<< "    Collisions                      = "
            << collisions << endl;
    }
    else
    {
        Info<< "    No collisions" << endl;
    }
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
