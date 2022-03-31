//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
/// \file B1PrimaryGeneratorAction.cc
/// \brief Implementation of the B1PrimaryGeneratorAction class

#include "B1PrimaryGeneratorAction.hh"

#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4RunManager.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4Event.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1PrimaryGeneratorAction::B1PrimaryGeneratorAction()
: G4VUserPrimaryGeneratorAction(),
  fParticleGun(nullptr)//,
  //fEnvelopeBox(0)
{
  G4int n_particle = 1;
  fParticleGun  = new G4ParticleGun(n_particle);

  //G4double envSizeXY = 0;
  //G4double envSizeZ = 0;

  //G4double size = 0.8;
  //G4double x0 = size * envSizeXY * (G4UniformRand()-0.5);
  //G4double y0 = size * envSizeXY * (G4UniformRand()-0.5);
  //G4double z0 = -0.5 * envSizeZ;

  G4double x0 = 0.;
  G4double y0 = 0.;
  G4double z0 = 0.;

  fParticleGun->SetParticlePosition(G4ThreeVector(x0,y0,z0));

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1PrimaryGeneratorAction::~B1PrimaryGeneratorAction()
{
  delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{

  //This function is called at the begining of ecah event
  //

  // In order to avoid dependence of PrimaryGeneratorAction
  // on DetectorConstruction class we get world volme
  // from G4LogicalVolumeStore.


  auto worldLV = G4LogicalVolumeStore::GetInstance()->GetVolume("World");

  // Check that the world volume has box shape
  G4Box* worldBox = nullptr;
  if (  worldLV ) {
    worldBox = dynamic_cast<G4Box*>(worldLV->GetSolid());
  }

  if ( worldBox ) {
    worldBox->GetZHalfLength();
  }
  else  {
    G4ExceptionDescription msg;
    msg << "World volume of box shape not found.\n";
    msg << "Perhaps you have changed geometry.\n";
    msg << "The gun will be place at the center.";
    G4Exception("B1PrimaryGeneratorAction::GeneratePrimaries()",
     "MyCode0002",JustWarning,msg);
  }

//------------------------------------------------





// Electron definitios as a particle

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();

  G4String particleName;
  G4ParticleDefinition* electron
    = particleTable->FindParticle(particleName="e-");
  fParticleGun->SetParticleDefinition(electron);


// Definition of the momentum direction of the particle emmited.
//Isotropic emission particle

  G4double Pi = 4. * atan(1.);
  G4double theta = 2 * Pi * G4UniformRand();
  G4double phi = Pi * G4UniformRand();

  G4double x_elec = sin(phi) * cos(theta);
  G4double y_elec = sin(phi) * sin(theta);
  G4double z_elec = cos(phi) ;

  G4ThreeVector pos_elec = G4ThreeVector( x_elec, y_elec, z_elec );

  fParticleGun->SetParticleMomentumDirection(pos_elec);

//Defiinition of the particle energy

  G4double Energy_elec = 4.0*MeV;
  fParticleGun->SetParticleEnergy(Energy_elec);

  fParticleGun->GeneratePrimaryVertex(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
