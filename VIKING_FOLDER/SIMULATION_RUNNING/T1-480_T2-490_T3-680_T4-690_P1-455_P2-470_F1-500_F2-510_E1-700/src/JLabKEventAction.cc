#include "JLabKEventAction.hh"

#include "JLabKRunAction.hh"
#include "JLabKVSteppingAction.hh"

#include "G4EventManager.hh"

JLabKEventAction::JLabKEventAction
(JLabKVSteppingAction* onePhotonSteppingAction)
: fpJLabKVSteppingAction(onePhotonSteppingAction)
{}

JLabKEventAction::~JLabKEventAction()
{}

void JLabKEventAction::BeginOfEventAction(const G4Event*)
{
  fpJLabKVSteppingAction->BeginOfEventAction();
}

void JLabKEventAction::EndOfEventAction(const G4Event*)
{
  fpJLabKVSteppingAction->EndOfEventAction();
}
