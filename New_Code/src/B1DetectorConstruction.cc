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
/// \file B1DetectorConstruction.cc
/// \brief Implementation of the B1DetectorConstruction class

#include "B1DetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4AssemblyVolume.hh"
#include "G4RotationMatrix.hh"
#include "G4PVReplica.hh"
#include "G4Transform3D.hh"



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1DetectorConstruction::B1DetectorConstruction()
: G4VUserDetectorConstruction(),
  fScoringVolume(0)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1DetectorConstruction::~B1DetectorConstruction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* B1DetectorConstruction::Construct()
{
// Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();


// Option to switch on/off checking of volumes overlaps
//
  G4bool checkOverlaps = true;


//******************************************************************************
// World
//******************************************************************************


  G4double world_sizeXY = 50.*cm;
  G4double world_sizeZ  = 50.*cm;
  G4Material* world_mat = nist->FindOrBuildMaterial("G4_Galactic");

  G4Box* solidWorld =
    new G4Box("World",                       //its name
       0.5*world_sizeXY, 0.5*world_sizeXY, 0.5*world_sizeZ);     //its size

  G4LogicalVolume* logicWorld =
    new G4LogicalVolume(solidWorld,          //its solid
                        world_mat,           //its material
                        "World");            //its name

  G4VPhysicalVolume* physWorld =
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(),       //at (0,0,0)
                      logicWorld,            //its logical volume
                      "World",               //its name
                      0,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      checkOverlaps);        //overlaps checking



//******************************************************************************
// Calorimeter
//******************************************************************************

  G4Material* cal_mat = nist->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");


  G4double cal_x = 10.*cm;
  G4double cal_y = 5.*cm;
  G4double cal_z = 5.*cm;

  G4Box* solidCal =
    new G4Box("Cal",
    cal_x/2, cal_y/2, cal_z/2);

  G4LogicalVolume* logicCal =
    new G4LogicalVolume(solidCal,         //its solid
                        cal_mat,          //its material
                        "Cal_LV");           //its name

//******************************************************************************
// Bars
//******************************************************************************


  G4Material* bar_mat = nist->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");


  G4double horBar_x = 0.2*cm;
  G4double horBar_y = cal_y;
  G4double horBar_z = cal_z/num_horBars;

  G4double verBar_x = 0.2*cm;
  G4double verBar_y = cal_y/num_verBars;
  G4double verBar_z = cal_z;


//------------------------------------------------------------------------------
// Horizontal bars
//------------------------------------------------------------------------------

  G4Box* solidHorBar =
    new G4Box("HorBar",
    horBar_x/2, horBar_y/2, horBar_z/2);

  G4LogicalVolume* logicHorBar =
    new G4LogicalVolume(solidHorBar,         //its solid
                        bar_mat,          //its material
                        "HorBar_LV");           //its name


//------------------------------------------------------------------------------
// Vertical bars
//------------------------------------------------------------------------------

  G4Box* solidVerBar =
    new G4Box("VerBar",
    verBar_x/2, verBar_y/2, verBar_z/2);

  G4LogicalVolume* logicVerBar =
    new G4LogicalVolume(solidVerBar,         //its solid
                        bar_mat,          //its material
                        "VerBar_LV");           //its name


 /*

//------------------------------------------------------------------------------
// Light guide
//------------------------------------------------------------------------------

  G4Material* lightGuide_mat = nist->FindOrBuildMaterial("G4_PLEXIGLASS");

  G4double lightGuide_x1 = 5.*cm;
  G4double lightGuide_x2 = 0.6*cm;
  G4double lightGuide_y1 = 5.*cm;
  G4double lightGuide_y2 = 0.6*cm;
  G4double lightGuide_z = 4.*cm;

  G4Trd* solidLightGuide =
  	new G4Trd("lightGuide",
      lightGuide_x1/2., 	//Half-length along X at the surface positioned at -dz
      lightGuide_x2/2.,	//Half-length along X at the surface positioned at +dz
      lightGuide_y1/2.,	//Half-length along Y at the surface positioned at -dz
      lightGuide_y2/2.,	//Half-length along Y at the surface positioned at +dz
      lightGuide_z/2.); //Half-length along Z axis

    G4LogicalVolume* logicLightGuide =
    new G4LogicalVolume(solidLightGuide,         //its solid
                        lightGuide_mat,          //its material
                        "lightGuide");           //its name

   G4ThreeVector pos_lightGuide = G4ThreeVector(dis_detsrc + bar_x + cal_x + lightGuide_z/2., 0.*cm, 0.*cm);

   G4RotationMatrix* yRot_lightGuide = new G4RotationMatrix; // Rotates Y and Z axes only
   yRot_lightGuide -> rotateY(-M_PI/2.*rad);	// Rotates 45 degrees

    new G4PVPlacement(yRot_lightGuide,                       //no rotation
                    pos_lightGuide,                    //at position
                    logicLightGuide,             //its logical volume
                    "lightGuide_LV",                //its name
                    logicWorld,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking
   */
   /*

  G4double dis_detsrc = 5.*cm;		//Source-Detector Distance

  G4ThreeVector pos_HorBar = G4ThreeVector(dis_detsrc + bar_x/2., 0.*cm, -( cal_z - bar_z )/2.);

  new G4PVPlacement(0,                       //no rotation
                    pos_HorBar,                    //at position
                    logicHorBar,             //its logical volume
                    "HorBar",                //its name
                    logicWorld,              //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking



 G4ThreeVector pos_Cal = G4ThreeVector(dis_detsrc + bar_x + cal_x/2., 0.*cm, 0*cm);

  new G4PVPlacement(0,                       //no rotation
                    pos_Cal,                    //at position
                    logicCal,             //its logical volume
                    "Cal",                //its name
                    logicWorld,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking


    */



//==============================================================================
//Assemblying the Clover
//==============================================================================


G4AssemblyVolume* assemblyCloverDetectors = new G4AssemblyVolume();

G4double detSouDis = 10.*cm;		//Source-Detector Distance

// Translation and Rotation of a clover inside the assembly

   G4ThreeVector T_Cal_0;
   G4ThreeVector T_Cal;
   G4RotationMatrix* R_Cal = new G4RotationMatrix;

   G4ThreeVector T_HorBar_0;
   G4ThreeVector T_HorBar;
   G4RotationMatrix* R_HorBar = new G4RotationMatrix;


   G4ThreeVector T_VerBar_0;
   G4ThreeVector T_VerBar;
   G4RotationMatrix* R_VerBar = new G4RotationMatrix;




 // Rotation of the assembly inside the world
    G4ThreeVector Tm;
    G4RotationMatrix* Rm = new G4RotationMatrix;


// Fill the assembly by the Clovers

   // Calorimeter
   T_Cal_0.setX( /*detSorDis*/ + horBar_x + verBar_x + cal_x/2. );
   T_Cal_0.setY( 0. );
   T_Cal_0.setZ( 0. );

   R_Cal->rotateX(0.*deg);
   R_Cal->rotateY(0.*deg);
   R_Cal->rotateY(0.*deg);

   //Horizontal bars
   T_HorBar_0.setX( /*detSorDis*/ + horBar_x/2. );
   T_HorBar_0.setY( 0. );
   T_HorBar_0.setZ( -( cal_z - horBar_z )/2. );

   //Vertical bars
   T_VerBar_0.setX( /*detSorDis*/ + horBar_x + verBar_x/2. );
   T_VerBar_0.setY( -( cal_y - verBar_y )/2. );
   T_VerBar_0.setZ( 0. );


G4double factor = 0.5 * (sqrt(num_telescopes) - 1.) ;
G4double theta[] {0., 72., 144., 216., 288};





for (int k = 0 ; k < sqrt(num_telescopes) ; k++ ){

	for(int j = 0; j < sqrt(num_telescopes) ; j++ ){


   		T_Cal.setX(T_Cal_0.getX());
   		T_Cal.setY(T_Cal_0.getY() + ((k-factor)*cal_y) );
   		T_Cal.setZ(T_Cal_0.getZ() + ((j-factor)*cal_z));

   		assemblyCloverDetectors->AddPlacedVolume( logicCal, T_Cal, R_Cal  );

   		for(int i=0 ; i < num_horBars ; i++){

   		T_HorBar.setX(T_HorBar_0.getX());
   		T_HorBar.setY(T_HorBar_0.getY() + ((k-factor)*cal_y));
   		T_HorBar.setZ(T_HorBar_0.getZ()+ i*horBar_z +  ((j-factor)*cal_z) );

   		R_HorBar->rotateX(0.*deg);
   		R_HorBar->rotateY(0.*deg);
   		R_HorBar->rotateZ(0.*deg);

   		assemblyCloverDetectors->AddPlacedVolume( logicHorBar, T_HorBar, R_HorBar  );
   		}

   		for(int i=0 ; i < num_verBars ; i++){
   		T_VerBar.setX(T_VerBar_0.getX());
   		T_VerBar.setY(T_VerBar_0.getY()+ i*verBar_y + ((k-factor)*cal_y));
   		T_VerBar.setZ(T_VerBar_0.getZ() +  ((j-factor)*cal_z));

   		R_VerBar->rotateX(0.*deg);
   		R_VerBar->rotateY(0.*deg);
   		R_VerBar->rotateZ(0.*deg);

   		assemblyCloverDetectors->AddPlacedVolume( logicVerBar, T_VerBar, R_VerBar  );
   		}

	}
}


    // Now instantiate the layers

   for( unsigned int i = 0; i < num_colvers ; i++ ) {
     // Translation of the assembly inside the world

     if(i == 0){
       Tm.setX(detSouDis*cos(theta[i]*deg));
       Tm.setY(detSouDis*sin(theta[i]*deg));
       Tm.setZ(0.*cm);

       Rm->rotateX(0.*deg);
       Rm->rotateY(0.*deg);
       Rm->rotateZ(theta[i]*deg);

     }else{
       Tm.setX(detSouDis*cos(theta[i]*deg));
       Tm.setY(detSouDis*sin(theta[i]*deg));
       Tm.setZ(0.*cm);

       Rm->rotateX(0.*deg);
       Rm->rotateY(0.*deg);
       Rm->rotateZ((theta[i] - theta[i-1])*deg);
     }





     assemblyCloverDetectors->MakeImprint( logicWorld, Tm, Rm );
   }


//============================================================================================






  // Set Calorimeter as scoring volume
  //
  fScoringVolume = logicCal;

  //
  //always return the physical World
  //
  return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
