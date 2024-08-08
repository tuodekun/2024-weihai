#include "PrimaryGeneratorAction.hh"

#include "G4RunManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include <vector>
#include "G4GenericMessenger.hh"

namespace B4
{

PrimaryGeneratorAction::PrimaryGeneratorAction()
// : fCurrentEnergy(1.*GeV),    // 默认能量
//    fNumParticlesPerEvent(1)   // 每个事件的默认粒子数
{
  G4int nofParticles = 1;
  fParticleGun = new G4ParticleGun(nofParticles);

  // Set up particle definition for a photon
  auto particleDefinition = G4ParticleTable::GetParticleTable()->FindParticle("gamma");
  fParticleGun->SetParticleDefinition(particleDefinition);

  // Array of energies [1GeV, 2GeV, 5GeV, 10GeV, 50GeV, 100GeV]
//  fEnergies = {1.*GeV, 2.*GeV, 5.*GeV, 10.*GeV, 50.*GeV, 100.*GeV};
//  fCurrentEnergyIndex = 0;
//  fPhotonCount = 0;
// Define commands via GenericMessenger
  fMessenger = new G4GenericMessenger(this, "/B4/generator/", "Primary generator control");
  fMessenger->DeclareProperty("energy", fCurrentEnergy, "Energy of primary particles");
  fMessenger->DeclareProperty("numParticlesPerEvent", fNumParticlesPerEvent, "Number of particles per event");
}

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fParticleGun;
  delete fMessenger;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  G4double worldZHalfLength = 0.;
  auto worldLV = G4LogicalVolumeStore::GetInstance()->GetVolume("World");

  G4Box* worldBox = nullptr;  
  if (worldLV) {
    worldBox = dynamic_cast<G4Box*>(worldLV->GetSolid());
  }

  if (worldBox) {
    worldZHalfLength = worldBox->GetZHalfLength();
  } else {
    G4ExceptionDescription msg;
    msg << "World volume of box shape not found." << G4endl;
    msg << "The gun will be placed in the center.";
    G4Exception("PrimaryGeneratorAction::GeneratePrimaries()",
      "MyCode0002", JustWarning, msg);
  }

  // Ensure the position is within the world volume but outside the detector
  G4double positionZ = -worldZHalfLength;
/*
  if (positionZ > worldZHalfLength) {
    positionZ = worldZHalfLength - 1.*mm; // Ensure it stays in the world volume
  }
*/

  // Generate random (theta, phi) within the specified ranges
  G4double theta = (2 + (4-2)*G4UniformRand()) * deg;
  G4double phi = (-CLHEP::pi + 2*CLHEP::pi*G4UniformRand());

  G4ThreeVector direction(std::sin(theta)*std::cos(phi), std::sin(theta)*std::sin(phi), std::cos(theta));

  fParticleGun->SetParticlePosition(G4ThreeVector(0., 0., positionZ));
  fParticleGun->SetParticleMomentumDirection(direction);

  fParticleGun->SetParticleEnergy(fCurrentEnergy);
  fParticleGun->GeneratePrimaryVertex(anEvent);

}
void PrimaryGeneratorAction::SetEnergy(G4double energy) {
  fCurrentEnergy = energy;
}

void PrimaryGeneratorAction::SetNumParticlesPerEvent(G4int num) {
  fNumParticlesPerEvent = num;
}

}

