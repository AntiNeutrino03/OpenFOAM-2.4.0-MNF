/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
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

Description

\*---------------------------------------------------------------------------*/

#include "dsmcReaction.H"
#include "IFstream.H"
#include "graph.H"
#include "dsmcCloud.H"


namespace Foam
{
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(dsmcReaction, 0);

defineRunTimeSelectionTable(dsmcReaction, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
dsmcReaction::dsmcReaction
(
    Time& t,
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    mesh_(cloud.mesh()),
    cloud_(cloud),
    nTotReactions_(0),
    reactWithLists_(false),
    propsDict_(dict.subDict(typeName + "Properties")),
    reactionName_(propsDict_.lookup("reactionName")),
    reactants_(propsDict_.lookup("reactants")),
    reactantIds_(),
    rDof1_(),
    rDof2_(),
    vDof1_(),
    vDof2_(),
    charge1_(),
    charge2_(),
    relax_(true),
    allowSplitting_(false),
    writeRatesToTerminal_(false),
    volume_(0.0),
    numberDensities_(2, 0.0)
{}

dsmcReaction::reactionProperties::reactionProperties
(
    dsmcCloud& cloud,
    dsmcParcel& p,
    dsmcParcel& q
)
:
    mesh_(cloud.mesh()),
    cloud_(cloud),
    typeIdP_(p.typeId()),
    typeIdQ_(q.typeId()),
    UP_(p.U()),
    UQ_(q.U()),
    ERotP_(p.ERot()),
    ERotQ_(q.ERot()),
    ELevelP_(p.ELevel()),
    ELevelQ_(q.ELevel()),
    vibLevelP_(p.vibLevel()),
    vibLevelQ_(q.vibLevel())
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

autoPtr<dsmcReaction> dsmcReaction::New
(
    Time& t,
    dsmcCloud& cloud,
    const dictionary& dict
)
{
    word dsmcReactionName
    (
        dict.lookup("reactionModel")
    );

    Info<< "Selecting the reaction model "
         << dsmcReactionName << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(dsmcReactionName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalError
            << "dsmcReaction::New(const dictionary&) : " << endl
            << "    unknown dsmc reaction model type "
            << dsmcReactionName
            << ", constructor not in hash table" << endl << endl
            << "    Valid reaction types are :" << endl;
        Info<< dictionaryConstructorTablePtr_->toc() << abort(FatalError);
    }

    return autoPtr<dsmcReaction>
    (
        cstrIter()(t, cloud, dict)
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dsmcReaction::~dsmcReaction()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void dsmcReaction::setCommonReactionProperties()
{
    // reading in reactants

//     reactants_ = propsDict_.lookup("reactants");

    if (reactants_.size() != 2)
    {
        FatalErrorIn("dsmcReaction::setCommonReactionProperties()")
            << "There should be two or more reactants, instead of "
            << reactants_.size() << nl
            << exit(FatalError);
    }
    
    reactantIds_.setSize(reactants_.size(), -1);

    allowSplitting_ = Switch(propsDict_.lookup("allowSplitting"));
    
    writeRatesToTerminal_ = Switch(propsDict_.lookup("writeRatesToTerminal"));

    forAll(reactantIds_, r)
    {
        reactantIds_[r] = findIndex(cloud_.typeIdList(), reactants_[r]);

        // check that reactants belong to the typeIdList
        if (reactantIds_[r] == -1)
        {
            FatalErrorIn("dsmcReaction::setCommonReactionProperties()")
                << "Cannot find type id: " << reactants_[r] << nl
                << exit(FatalError);
        }
    }
    
    rDof1_ =
        cloud_.constProps(reactantIds_[0]).rotationalDegreesOfFreedom();
        
    rDof2_ =
        cloud_.constProps(reactantIds_[1]).rotationalDegreesOfFreedom();
        
    vDof1_ =
        cloud_.constProps(reactantIds_[0]).vibrationalDegreesOfFreedom();
        
    vDof2_ =
        cloud_.constProps(reactantIds_[1]).vibrationalDegreesOfFreedom();
        
    charge1_ = cloud_.constProps(reactantIds_[0]).charge();
    
    charge2_ = cloud_.constProps(reactantIds_[1]).charge();

}

label dsmcReaction::nTotReactions() const
{
    return nTotReactions_;
}


label dsmcReaction::nReactionsPerTimeStep() const
{
    return nReactionsPerTimeStep_;
}


label& dsmcReaction::nReactionsPerTimeStep()
{
    return nReactionsPerTimeStep_;
}


const dsmcCloud& dsmcReaction::cloud() const
{
    return cloud_;
}


bool dsmcReaction::reactWithLists() const
{
    return reactWithLists_;
}

vector& dsmcReaction::reactionProperties::UP()
{
    return  UP_;
}

vector& dsmcReaction::reactionProperties::UQ()
{
    return  UQ_;
}

scalar& dsmcReaction::reactionProperties::ERotP()
{
    return ERotP_;
}

scalar& dsmcReaction::reactionProperties::ERotQ()
{
    return ERotQ_;
}

scalar dsmcReaction::reactionProperties::EVibP() const
{
    if(cloud_.constProps(typeIdP_).vibrationalDegreesOfFreedom() > VSMALL)
    {
        return vibLevelP_[0]*cloud_.constProps(typeIdP_).thetaV()[0]
                         *physicoChemical::k.value();
    }
    else
    {
        return 0.0;
    }
}

scalar dsmcReaction::reactionProperties::EVibQ() const
{
    if(cloud_.constProps(typeIdQ_).vibrationalDegreesOfFreedom() > VSMALL)
    {
        return vibLevelQ_[0]*cloud_.constProps(typeIdQ_).thetaV()[0]
                         *physicoChemical::k.value();
    }
    else
    {
        return 0.0;
    }
}

scalar dsmcReaction::reactionProperties::EEleP() const
{
    return cloud_.constProps(typeIdP_).electronicEnergyList()[ELevelP_];
}

scalar dsmcReaction::reactionProperties::EEleQ() const
{
    return cloud_.constProps(typeIdQ_).electronicEnergyList()[ELevelQ_];
}

scalar dsmcReaction::reactionProperties::mP() const
{
    return cloud_.constProps(typeIdP_).mass();
}

scalar dsmcReaction::reactionProperties::mQ() const
{
    return  cloud_.constProps(typeIdQ_).mass();
}

scalar dsmcReaction::reactionProperties::mR() const
{
    return  cloud_.constProps(typeIdP_).mass()*cloud_.constProps(typeIdQ_).mass()/(cloud_.constProps(typeIdP_).mass() + cloud_.constProps(typeIdQ_).mass());
}

scalar dsmcReaction::reactionProperties::translationalEnergy() const
{
    scalar translationalEnergy = 0.5*(cloud_.constProps(typeIdP_).mass()* cloud_.constProps(typeIdQ_).mass()/(cloud_.constProps(typeIdP_).mass() +  cloud_.constProps(typeIdQ_).mass()))*magSqr(UP_ - UQ_);
    
    return  translationalEnergy;
}

scalar dsmcReaction::reactionProperties::omegaPQ() const
{
    return  0.5*(cloud_.constProps(typeIdP_).omega() + cloud_.constProps(typeIdQ_).omega());
}

scalar dsmcReaction::reactionProperties::thetaDP() const
{
    return  cloud_.constProps(typeIdP_).thetaD()[0];
}

scalar dsmcReaction::reactionProperties::thetaDQ() const
{
    return  cloud_.constProps(typeIdQ_).thetaD()[0];
}

scalar dsmcReaction::reactionProperties::thetaVP() const
{
    return  cloud_.constProps(typeIdP_).thetaV()[0];
}

scalar dsmcReaction::reactionProperties::thetaVQ() const
{
    return  cloud_.constProps(typeIdQ_).thetaV()[0];
}

scalar dsmcReaction::reactionProperties::ZrefP() const
{
    return  cloud_.constProps(typeIdP_).Zref()[0];
}

scalar dsmcReaction::reactionProperties::ZrefQ() const
{
    return  cloud_.constProps(typeIdQ_).Zref()[0];
}

scalar dsmcReaction::reactionProperties::refTempZvP() const
{
    return  cloud_.constProps(typeIdP_).TrefZv()[0];
}

scalar dsmcReaction::reactionProperties::refTempZvQ() const
{
    return  cloud_.constProps(typeIdQ_).TrefZv()[0];
}

label dsmcReaction::reactionProperties::charDissLevelP() const
{
    return  cloud_.constProps(typeIdP_).charDissQuantumLevel()[0];
}

label dsmcReaction::reactionProperties::charDissLevelQ() const
{
    return  cloud_.constProps(typeIdQ_).charDissQuantumLevel()[0];
}

label dsmcReaction::reactionProperties::jMaxP() const
{
    return  cloud_.constProps(typeIdP_).numberOfElectronicLevels()-1;
}

label dsmcReaction::reactionProperties::jMaxQ() const
{
    return  cloud_.constProps(typeIdQ_).numberOfElectronicLevels()-1;
}


label dsmcReaction::reactionProperties::rotationalDofP() const
{
    return   cloud_.constProps(typeIdP_).rotationalDegreesOfFreedom();
}

label dsmcReaction::reactionProperties::rotationalDofQ() const
{
    return   cloud_.constProps(typeIdQ_).rotationalDegreesOfFreedom();
}

vector dsmcReaction::reactionProperties::Ucm() const
{
    scalar mP = cloud_.constProps(typeIdP_).mass();
    scalar mQ = cloud_.constProps(typeIdQ_).mass();
    
    return (mP*UP_ + mQ*UQ_)/(mP + mQ);
}

const List<scalar>& dsmcReaction::reactionProperties::EEListP() const
{
    return  cloud_.constProps(typeIdP_).electronicEnergyList();
}

const List<scalar>& dsmcReaction::reactionProperties::EEListQ() const
{
    return  cloud_.constProps(typeIdQ_).electronicEnergyList();
}

const List<label>& dsmcReaction::reactionProperties::gListP() const
{
    return  cloud_.constProps(typeIdP_).degeneracyList();
}

const List<label>& dsmcReaction::reactionProperties::gListQ() const
{
    return  cloud_.constProps(typeIdQ_).degeneracyList();
}

void dsmcReaction::reactionProperties::dissociateP
(
    const scalar& heatOfReaction,
    const List<label>& productIds,
    dsmcParcel& p,
    dsmcParcel& q
)
{
    scalar ChiB = 2.5 - this->omegaPQ();

    scalar heatOfReactionDissociationJoules =
                            heatOfReaction*physicoChemical::k.value();

    scalar translationalEnergy = this->translationalEnergy();

    label ELevelQ = 0;

    translationalEnergy += heatOfReactionDissociationJoules + this->EVibP();

    translationalEnergy += this->EEleQ();

    ELevelQ = cloud_.postCollisionElectronicEnergyLevel
                    (
                        translationalEnergy,
                        (this->jMaxQ() + 1),
                        this->omegaPQ(),
                        this->EEListQ(),
                        this->gListQ()
                    );

    translationalEnergy -= this->EEListQ()[ELevelQ];

    label vibLevelQ = 0;
    scalar ERotQ = 0.0;

    if (this->rotationalDofQ() > VSMALL)
    {
        translationalEnergy += this->EVibQ();

        label iMax = translationalEnergy
                        / (physicoChemical::k.value()*this->thetaVQ());

        vibLevelQ = cloud_.postCollisionVibrationalEnergyLevel
                        (
                                true,
                                q.vibLevel()[0],
                                iMax,
                                this->thetaVQ(),
                                this->thetaDQ(),
                                this->refTempZvQ(),
                                this->omegaPQ(),
                                this->ZrefQ(),
                                translationalEnergy
                            );

        translationalEnergy -= vibLevelQ
                            *this->thetaVQ()*physicoChemical::k.value();

        translationalEnergy += this->ERotQ();

        ERotQ = translationalEnergy
        *cloud_.postCollisionRotationalEnergy(this->rotationalDofQ(), ChiB);

        translationalEnergy -= ERotQ;
    }

    scalar relVel = sqrt(2.0*translationalEnergy/this->mR());

    // Variable Hard Sphere collision part

    scalar cosTheta = 2.0*cloud_.rndGen().scalar01() - 1.0;

    scalar sinTheta = sqrt(1.0 - cosTheta*cosTheta);

    scalar phi = twoPi*cloud_.rndGen().scalar01();

    vector postCollisionRelU =
        relVel
        *vector
            (
                cosTheta,
                sinTheta*cos(phi),
                sinTheta*sin(phi)
            );

    vector UP = this->Ucm() + (postCollisionRelU*this->mQ()/(this->mP() + this->mQ()));
    vector UQ = this->Ucm() - (postCollisionRelU*this->mP()/(this->mP() + this->mQ()));

    label typeId1 = productIds[0];
    label typeId2 = productIds[1];

    // Mass of Product one and two
    const scalar& mP1 = cloud_.constProps(typeId1).mass();
    const scalar& mP2 = cloud_.constProps(typeId2).mass();

    scalar mRatoms = mP1*mP2/(mP1 + mP2);

    translationalEnergy = this->ERotP() + this->EEleP();

    scalar cRatoms = sqrt(2.0*translationalEnergy/mRatoms);

    // Variable Hard Sphere collision part
    scalar cosTheta2 = 2.0*cloud_.rndGen().scalar01() - 1.0;

    scalar sinTheta2 = sqrt(1.0 - cosTheta2*cosTheta2);

    scalar phi2 = twoPi*cloud_.rndGen().scalar01();

    vector postCollisionRelU2 = cRatoms
    *vector
        (
            cosTheta2,
            sinTheta2*cos(phi2),
            sinTheta2*sin(phi2)
        );

    vector UP1 = UP + (postCollisionRelU2)*mP2/(mP1 + mP2);
    vector UP2 = UP - (postCollisionRelU2)*mP1/(mP1 + mP2);

    q.U() = UQ;
    q.ERot() = ERotQ;
    if (this->rotationalDofQ() > VSMALL)
    {
        q.vibLevel()[0] = vibLevelQ;
    }

    q.ELevel() = ELevelQ;

    // Molecule P will dissociate
    const vector& position(p.position());

    label cell = -1;
    label tetFace = -1;
    label tetPt = -1;

    mesh_.findCellFacePt
    (
        position,
        cell,
        tetFace,
        tetPt
    );

    scalar RWF = p.RWF();
    label classification = p.classification();
    scalarField wallTemperature(3, 0.0);
    vectorField wallVectors(3, vector::zero);
    labelList vibLevel(0, 0);

    p.position() = position;
    p.U() = UP1;
    p.RWF() = RWF;
    p.ERot() = 0.0;
    p.ELevel() = 0;
    p.cell() = cell;
    p.tetFace() = tetFace;
    p.tetPt() = tetPt;
    p.typeId() = typeId1;
    p.newParcel() = 0;
    p.classification() = classification;
    p.stuckToWall() = 0;
    p.wallTemperature() = wallTemperature;
    p.wallVectors() = wallVectors;
    p.vibLevel()= vibLevel;
    
    // insert new product
    cloud_.addNewParcel
    (
        position,
        UP2,
        position,
        RWF,
        0.0,
        0,
        cell,
        tetFace,
        tetPt,
        typeId2,
        0,
        classification,
        0,
        wallTemperature,
        wallVectors,
        vibLevel
    );
}

void dsmcReaction::reactionProperties::dissociateQ
(
    const scalar& heatOfReaction,
    const List<label>& productIds,
    dsmcParcel& p,
    dsmcParcel& q
)
{
    scalar ChiB = 2.5 - this->omegaPQ();
    
    scalar heatOfReactionDissociationJoules =
            heatOfReaction*physicoChemical::k.value();
            
    scalar translationalEnergy = this->translationalEnergy();

    translationalEnergy += heatOfReactionDissociationJoules + this->EVibQ();

    translationalEnergy += this->EEleP();

    label ELevelP = cloud_.postCollisionElectronicEnergyLevel
                (
                    translationalEnergy,
                    (this->jMaxP() + 1),
                    this->omegaPQ(),
                    this->EEListP(),
                    this->gListP()
                );

    translationalEnergy -= this->EEListP()[ELevelP];

    label vibLevelP = 0;
    scalar ERotP = 0.0;

    if (this->rotationalDofP() > VSMALL)
    {
        translationalEnergy += this->EVibP();

        label iMax = translationalEnergy
                        / (physicoChemical::k.value()*this->thetaVP());

        vibLevelP = cloud_.postCollisionVibrationalEnergyLevel
                        (
                                true,
                                p.vibLevel()[0],
                                iMax,
                                this->thetaVP(),
                                this->thetaDP(),
                                this->refTempZvP(),
                                this->omegaPQ(),
                                this->ZrefP(),
                                translationalEnergy
                            );

        translationalEnergy -= vibLevelP*this->thetaVP()
                                *physicoChemical::k.value();

        translationalEnergy += this->ERotP();

        ERotP = translationalEnergy
        *cloud_.postCollisionRotationalEnergy(this->rotationalDofP(), ChiB);

        translationalEnergy -= ERotP;
    }

    scalar relVelNonDissoMol = sqrt(2.0*translationalEnergy/this->mR());

    // Variable Hard Sphere collision part

    scalar cosTheta = 2.0*cloud_.rndGen().scalar01() - 1.0;

    scalar sinTheta = sqrt(1.0 - cosTheta*cosTheta);

    scalar phi = twoPi*cloud_.rndGen().scalar01();

    vector postCollisionRelU =
    relVelNonDissoMol
    *vector
        (
            cosTheta,
            sinTheta*cos(phi),
            sinTheta*sin(phi)
        );

    vector UP = this->Ucm() + (postCollisionRelU*this->mQ()/(this->mP() + this->mQ()));
    vector UQ = this->Ucm() - (postCollisionRelU*this->mP()/(this->mP() + this->mQ()));

    const label typeId1 = productIds[0];
    const label typeId2 = productIds[1];

    // Mass of Product one and two
    scalar mP1 = cloud_.constProps(typeId1).mass();
    scalar mP2 = cloud_.constProps(typeId2).mass();

    scalar mRatoms = mP1*mP2/(mP1 + mP2);

    translationalEnergy = this->ERotQ() + this->EEleQ();

    scalar cRatoms = sqrt(2.0*translationalEnergy/mRatoms);

    // Variable Hard Sphere collision part
    scalar cosTheta2 = 2.0*cloud_.rndGen().scalar01() - 1.0;

    scalar sinTheta2 = sqrt(1.0 - cosTheta2*cosTheta2);

    scalar phi2 = twoPi*cloud_.rndGen().scalar01();

    vector postCollisionRelU2 = cRatoms
    *vector
    (
        cosTheta2,
        sinTheta2*cos(phi2),
        sinTheta2*sin(phi2)
    );


    vector uQ1 = UQ + postCollisionRelU2*mP2/(mP1 + mP2);
    vector uQ2 = UQ - postCollisionRelU2*mP1/(mP1 + mP2);

    p.U() = UP;
    p.ERot() = ERotP;
    if (this->rotationalDofP() > VSMALL)
    {
        p.vibLevel()[0] = vibLevelP;
    }
    p.ELevel() = ELevelP;
    
    // Molecule Q will dissociate
    const vector& position(q.position());

    label cell = -1;
    label tetFace = -1;
    label tetPt = -1;

    mesh_.findCellFacePt
    (
        position,
        cell,
        tetFace,
        tetPt
    );

    scalar RWF = q.RWF();
    label classification = q.classification();
    scalarField wallTemperature(3, 0.0);
    vectorField wallVectors(3, vector::zero);
    labelList vibLevel(0, 0);

    q.position() = position;
    q.U() = uQ1;
    q.RWF() = RWF;
    q.ERot() = 0.0;
    q.ELevel() = 0;
    q.cell() = cell;
    q.tetFace() = tetFace;
    q.tetPt() = tetPt;
    q.typeId() = typeId1;
    q.newParcel() = 0;
    q.classification() = classification;
    q.stuckToWall() = 0;
    q.wallTemperature() = wallTemperature;
    q.wallVectors() = wallVectors;
    q.vibLevel() = vibLevel;
    
    // insert new product
    cloud_.addNewParcel
    (
        position,
        uQ2,
        position,
        RWF,
        0.0,
        0,
        cell,
        tetFace,
        tetPt,
        typeId2,
        0,
        classification,
        0,
        wallTemperature,
        wallVectors,
        vibLevel
    );
}

void dsmcReaction::reactionProperties::ioniseP
(
    const scalar& heatOfReaction,
    const List<label>& productIds,
    dsmcParcel& p,
    dsmcParcel& q
)
{
    scalar ChiB = 2.5 - this->omegaPQ();
    
    scalar heatOfReactionIonisationJoules = heatOfReaction*physicoChemical::k.value();
                
    scalar translationalEnergy = this->translationalEnergy();

    translationalEnergy += heatOfReactionIonisationJoules + this->EEleP();

    translationalEnergy += this->EEleQ();

    label ELevelQ = cloud_.postCollisionElectronicEnergyLevel
                    (
                        translationalEnergy,
                        (this->jMaxQ() + 1),
                        this->omegaPQ(),
                        this->EEListQ(),
                        this->gListQ()
                    );

    translationalEnergy -= this->EEListQ()[ELevelQ];

    label vibLevelQ = 0;
    scalar ERotQ = 0.0;

    if (this->rotationalDofQ() > VSMALL)
    {
        translationalEnergy += this->EVibQ();

        label iMax = translationalEnergy
                        / (physicoChemical::k.value()*this->thetaVQ());

        vibLevelQ = cloud_.postCollisionVibrationalEnergyLevel
                        (
                                true,
                                q.vibLevel()[0],
                                iMax,
                                this->thetaVQ(),
                                this->thetaDQ(),
                                this->refTempZvQ(),
                                this->omegaPQ(),
                                this->ZrefQ(),
                                translationalEnergy
                            );

        translationalEnergy -= vibLevelQ
                        *this->thetaVQ()*physicoChemical::k.value();

        translationalEnergy += this->ERotQ();

        ERotQ = translationalEnergy
        *cloud_.postCollisionRotationalEnergy(this->rotationalDofQ(), ChiB);

        translationalEnergy -= ERotQ;
    }

    scalar relVelNonDissoMol = sqrt(2.0*translationalEnergy/this->mR());

    // Variable Hard Sphere collision part

    scalar cosTheta = 2.0*cloud_.rndGen().scalar01() - 1.0;

    scalar sinTheta = sqrt(1.0 - cosTheta*cosTheta);

    scalar phi = twoPi*cloud_.rndGen().scalar01();

    vector postCollisionRelU =
        relVelNonDissoMol
        *vector
            (
                cosTheta,
                sinTheta*cos(phi),
                sinTheta*sin(phi)
            );

    vector UP = this->Ucm() + (postCollisionRelU*this->mQ()/(this->mP() + this->mQ()));
    vector UQ = this->Ucm() - (postCollisionRelU*this->mP()/(this->mP() + this->mQ()));

    const label typeId1 = productIds[0];
    const label typeId2 = productIds[1];

    // Mass of Product one and two
    scalar mP1 = cloud_.constProps(typeId1).mass();
    scalar mP2 = cloud_.constProps(typeId2).mass();

    scalar mRatoms = mP1*mP2/(mP1 + mP2);

    translationalEnergy = this->ERotP() + this->EVibP();

    scalar cRatoms = sqrt(2.0*translationalEnergy/mRatoms);

    // Variable Hard Sphere collision part
    scalar cosTheta2 = 2.0*cloud_.rndGen().scalar01() - 1.0;

    scalar sinTheta2 = sqrt(1.0 - cosTheta2*cosTheta2);

    scalar phi2 = twoPi*cloud_.rndGen().scalar01();

    vector postCollisionRelU2 = cRatoms
    *vector
        (
            cosTheta2,
            sinTheta2*cos(phi2),
            sinTheta2*sin(phi2)
        );


    vector uP1 = UP + postCollisionRelU2*mP2/(mP1 + mP2);
    vector uP2 = UP - postCollisionRelU2*mP1/(mP1 + mP2);

    q.U() = UQ;
    q.ERot() = ERotQ;
    if (this->rotationalDofQ() > VSMALL)
    {
        q.vibLevel()[0] = vibLevelQ;
    }
    q.ELevel() = ELevelQ;

    vector position(p.position());

    label cell = -1;
    label tetFace = -1;
    label tetPt = -1;

    mesh_.findCellFacePt
    (
        position,
        cell,
        tetFace,
        tetPt
    );

    scalar RWF = p.RWF();
    label classification = p.classification();
    scalarField wallTemperature(3, 0.0);
    vectorField wallVectors(3, vector::zero);
    labelList vibLevel(0, 0);
    labelList electronVibLevel(1, 0);
    electronVibLevel[0] = 0;
    
    if(this->rotationalDofP() > VSMALL)
    {
        vibLevel.resize(1);
        vibLevel[0] = 0;
    }
    
    p.position() = position;
    p.U() = uP1;
    p.RWF() = RWF;
    p.ERot() = 0.0;
    p.ELevel() = 0;
    p.cell() = cell;
    p.tetFace() = tetFace;
    p.tetPt() = tetPt;
    p.typeId() = typeId1;
    p.newParcel() = 0;
    p.classification() = classification;
    p.stuckToWall() = 0;
    p.wallTemperature() = wallTemperature;
    p.wallVectors() = wallVectors;
    p.vibLevel()= vibLevel;
    
    // insert new product (electron)
    cloud_.addNewParcel
    (
        position,
        uP2,
        position,
        RWF,
        0.0,
        0,
        cell,
        tetFace,
        tetPt,
        typeId2,
        0,
        classification,
        0,
        wallTemperature,
        wallVectors,
        electronVibLevel
    );
}

void dsmcReaction::reactionProperties::ioniseQ
(
    const scalar& heatOfReaction,
    const List<label>& productIds,
    dsmcParcel& p,
    dsmcParcel& q
)
{
    scalar ChiB = 2.5 - this->omegaPQ();
    
    scalar heatOfReactionIonisationJoules =
            heatOfReaction*physicoChemical::k.value();
                            
    scalar translationalEnergy = this->translationalEnergy();

    translationalEnergy += heatOfReactionIonisationJoules + this->EEleQ();

    translationalEnergy += this->EEleP();

    label ELevelP = cloud_.postCollisionElectronicEnergyLevel
                    (
                        translationalEnergy,
                        (this->jMaxP() + 1),
                        this->omegaPQ(),
                        this->EEListP(),
                        this->gListP()
                    );

    translationalEnergy -= this->EEListP()[ELevelP];

    label vibLevelP = 0;
    scalar ERotP = 0.0;

    if (this->rotationalDofP() > VSMALL)
    {
        translationalEnergy += this->EVibP();

        label iMax = (translationalEnergy
                        / (physicoChemical::k.value()*this->thetaVP()));

        vibLevelP = cloud_.postCollisionVibrationalEnergyLevel
                        (
                                true,
                                p.vibLevel()[0],
                                iMax,
                                this->thetaVP(),
                                this->thetaDP(),
                                this->refTempZvP(),
                                this->omegaPQ(),
                                this->ZrefP(),
                                translationalEnergy
                            );

        translationalEnergy -= vibLevelP*this->thetaVP()
                            *physicoChemical::k.value();

        translationalEnergy += this->ERotP();

        ERotP = translationalEnergy
        *cloud_.postCollisionRotationalEnergy(this->rotationalDofP(), ChiB);

        translationalEnergy -= ERotP;
    }

    scalar relVelNonDissoMol = sqrt(2.0*translationalEnergy/this->mR());

    // Variable Hard Sphere collision part

    scalar cosTheta = 2.0*cloud_.rndGen().scalar01() - 1.0;

    scalar sinTheta = sqrt(1.0 - cosTheta*cosTheta);

    scalar phi = twoPi*cloud_.rndGen().scalar01();

    vector postCollisionRelU =
        relVelNonDissoMol
        *vector
            (
                cosTheta,
                sinTheta*cos(phi),
                sinTheta*sin(phi)
            );

    vector UP = this->Ucm() + (postCollisionRelU*this->mQ()/(this->mP() + this->mQ()));
    vector UQ = this->Ucm() - (postCollisionRelU*this->mP()/(this->mP() + this->mQ()));

    const label typeId1 = productIds[0];
    const label typeId2 = productIds[1];

    // Mass of Product one and two
    scalar mP1 = cloud_.constProps(typeId1).mass();
    scalar mP2 = cloud_.constProps(typeId2).mass();

    scalar mRatoms = mP1*mP2/(mP1 + mP2);

    translationalEnergy = this->ERotQ() + this->EVibQ();

    scalar cRatoms = sqrt(2.0*translationalEnergy/mRatoms);

    // Variable Hard Sphere collision part
    scalar cosTheta2 = 2.0*cloud_.rndGen().scalar01() - 1.0;

    scalar sinTheta2 = sqrt(1.0 - cosTheta2*cosTheta2);

    scalar phi2 = twoPi*cloud_.rndGen().scalar01();

    vector postCollisionRelU2 = cRatoms
    *vector
        (
            cosTheta2,
            sinTheta2*cos(phi2),
            sinTheta2*sin(phi2)
        );


    vector uQ1 = UQ + postCollisionRelU2*mP2/(mP1 + mP2);
    vector uQ2 = UQ - postCollisionRelU2*mP1/(mP1 + mP2);

    p.U() = UP;
    p.ERot() = ERotP;
    if (this->rotationalDofP() > VSMALL)
    {
        p.vibLevel()[0] = vibLevelP;
    }
    p.ELevel() = ELevelP;

    // Molecule Q will ionise.
    vector position(q.position());

    label cell = -1;
    label tetFace = -1;
    label tetPt = -1;

    mesh_.findCellFacePt
    (
        position,
        cell,
        tetFace,
        tetPt
    );

    scalar RWF = q.RWF();
    label classification = q.classification();
    scalarField wallTemperature(3, 0.0);
    vectorField wallVectors(3, vector::zero);
    labelList vibLevel(0, 0);
    labelList electronVibLevel(1, 0);
    electronVibLevel[0] = 0;
    
    if(this->rotationalDofQ() > VSMALL)
    {
        vibLevel.resize(1);
        vibLevel[0] = 0;
    }
    
    q.position() = position;
    q.U() = uQ1;
    q.RWF() = RWF;
    q.ERot() = 0.0;
    q.ELevel() = 0;
    q.cell() = cell;
    q.tetFace() = tetFace;
    q.tetPt() = tetPt;
    q.typeId() = typeId1;
    q.newParcel() = 0;
    q.classification() = classification;
    q.stuckToWall() = 0;
    q.wallTemperature() = wallTemperature;
    q.wallVectors() = wallVectors;
    q.vibLevel() = vibLevel;

    // insert new product (electron)
    cloud_.addNewParcel
    (
        position,
        uQ2,
        position,
        RWF,
        0.0,
        0,
        cell,
        tetFace,
        tetPt,
        typeId2,
        0,
        classification,
        0,
        wallTemperature,
        wallVectors,
        electronVibLevel
    );
}

void dsmcReaction::reactionProperties::associativeIonisation
(
    const scalar& heatOfReactionIntermediateIonisation,
    const scalar& heatOfReactionRecombination,
    const List<label>& assIonProductIds,
    dsmcParcel& p,
    dsmcParcel& q
)
{
    scalar heatOfReactionIonisationJoules =
                    heatOfReactionIntermediateIonisation
                    *physicoChemical::k.value();
    scalar heatOfReactionRecombinationJoules =
                    heatOfReactionRecombination
                    *physicoChemical::k.value();

    scalar translationalEnergy = this->translationalEnergy();

    translationalEnergy += heatOfReactionRecombinationJoules
                            + heatOfReactionIonisationJoules;

    translationalEnergy += this->EEleP() + this->EEleQ() + this->ERotP() 
                         + this->ERotQ() + this->EVibP() + this->EVibQ();

    // centre of mass velocity
    vector Ucm = (this->mP()*this->UP() + this->mQ()*this->UQ())/
                                                    (this->mP() + this->mQ());
      
    label ELevelNewP = 0;
    label ELevelNewQ = 0;
    label vibLevelNewP = 0;
    scalar ERotNewP = 0.0;
    scalar ERotNewQ = 0.0;

    scalar jMaxNewP =
        cloud_.constProps(assIonProductIds[0]).numberOfElectronicLevels();
    scalar jMaxNewQ =
        cloud_.constProps(assIonProductIds[1]).numberOfElectronicLevels();
    
    const scalarList& EElistNewP =
        cloud_.constProps(assIonProductIds[0]).electronicEnergyList();
    const scalarList& EElistNewQ =
        cloud_.constProps(assIonProductIds[1]).electronicEnergyList();
        
    const labelList& gListNewP =
        cloud_.constProps(assIonProductIds[0]).degeneracyList();
    const labelList& gListNewQ =
        cloud_.constProps(assIonProductIds[1]).degeneracyList();

    scalar omegaNewPQ =
        0.5
        *(
            cloud_.constProps(assIonProductIds[0]).omega()
            + cloud_.constProps(assIonProductIds[1]).omega()
        );
    
    ELevelNewP = cloud_.postCollisionElectronicEnergyLevel
                    (
                        translationalEnergy,
                        jMaxNewP,
                        omegaNewPQ,
                        EElistNewP,
                        gListNewP
                    );

    translationalEnergy -= EElistNewP[ELevelNewP];
    
    ELevelNewQ = cloud_.postCollisionElectronicEnergyLevel
                    (
                        translationalEnergy,
                        jMaxNewQ,
                        omegaNewPQ,
                        EElistNewQ,
                        gListNewQ
                    );

    translationalEnergy -= EElistNewQ[ELevelNewQ];
    
    scalar rotationalDofNewP =
            cloud_.constProps(assIonProductIds[0]).rotationalDegreesOfFreedom();

    if(rotationalDofNewP > VSMALL)
    {
        scalar thetaVNewP =
                cloud_.constProps(assIonProductIds[0]).thetaV()[0];
        scalar thetaDNewP =
                cloud_.constProps(assIonProductIds[0]).thetaD()[0];
                
        scalar ZrefNewP = cloud_.constProps(assIonProductIds[0]).Zref()[0];
        scalar refTempZvNewP =
                cloud_.constProps(assIonProductIds[0]).TrefZv()[0];
        
        scalar ChiB = 2.5 - omegaNewPQ;


        if (rotationalDofNewP > VSMALL)
        {
            label iMax = translationalEnergy
                        / (physicoChemical::k.value()*thetaVNewP);

            vibLevelNewP = cloud_.postCollisionVibrationalEnergyLevel
                            (
                                    true,
                                    0,
                                    iMax,
                                    thetaVNewP,
                                    thetaDNewP,
                                    refTempZvNewP,
                                    omegaNewPQ,
                                    ZrefNewP,
                                    translationalEnergy
                                );

            translationalEnergy -= vibLevelNewP*thetaVNewP
                                    *physicoChemical::k.value();

            ERotNewP = translationalEnergy
                        *cloud_.postCollisionRotationalEnergy(
                            rotationalDofNewP, ChiB);

            translationalEnergy -= ERotNewP;
        }
    }

    scalar mP2 = cloud_.constProps(assIonProductIds[0]).mass();
    scalar mQ2 = cloud_.constProps(assIonProductIds[1]).mass();

    scalar mR = mP2*mQ2/(mP2 + mQ2);

    scalar relVel = sqrt((2.0*translationalEnergy)/mR);

    // Variable Hard Sphere collision part for collision
    scalar cosTheta = 2.0*cloud_.rndGen().scalar01() - 1.0;

    scalar sinTheta = sqrt(1.0 - cosTheta*cosTheta);

    scalar phi = twoPi*cloud_.rndGen().scalar01();

    vector postCollisionRelU =
        relVel
        *vector
            (
                cosTheta,
                sinTheta*cos(phi),
                sinTheta*sin(phi)
            );

    vector UP = Ucm + (postCollisionRelU*mQ2/(mP2 + mQ2));
    vector UQ = Ucm - (postCollisionRelU*mP2/(mP2 + mQ2));
    
    vector positionP(p.position());
    
    label cell = -1;
    label tetFace = -1;
    label tetPt = -1;

    mesh_.findCellFacePt
    (
        positionP,
        cell,
        tetFace,
        tetPt
    );

    scalar RWFp = p.RWF();
    labelList vibLevelP(0, 0);
    if(rotationalDofNewP > VSMALL)
    {
        vibLevelP.resize(1);
        vibLevelP[0] = vibLevelNewP;
    }

    p.position() = positionP;
    p.U() = UP;
    p.RWF() = RWFp;
    p.ERot() = ERotNewP;
    p.ELevel() = ELevelNewP;
    p.cell() = cell;
    p.tetFace() = tetFace;
    p.tetPt() = tetPt;
    p.typeId() = assIonProductIds[0];
    p.newParcel() = 0;
    p.vibLevel() = vibLevelP;
    
    vector positionQ(q.position());
    
    scalar RWFq = q.RWF();
    labelList vibLevelQ(0, 0);
        
    cell = -1;
    tetFace = -1;
    tetPt = -1;

    mesh_.findCellFacePt
    (
        positionQ,
        cell,
        tetFace,
        tetPt
    );
    
    q.position() = positionQ;
    q.U() = UQ;
    q.RWF() = RWFq;
    q.ERot() = ERotNewQ;
    q.ELevel() = ELevelNewQ;
    q.cell() = cell;
    q.tetFace() = tetFace;
    q.tetPt() = tetPt;
    q.typeId() = assIonProductIds[1];
    q.newParcel() = 0;
    q.vibLevel() = vibLevelQ;
}

} // End namespace Foam

// ************************************************************************* //
