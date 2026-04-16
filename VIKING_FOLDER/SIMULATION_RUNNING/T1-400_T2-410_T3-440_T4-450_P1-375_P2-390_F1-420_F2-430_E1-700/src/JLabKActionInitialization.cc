#include "JLabKActionInitialization.hh"

#include "JLabKPrimaryGeneratorAction.hh"
#include "JLabKRunAction.hh"
#include "JLabKEventAction.hh"
#include "JLabKSteppingAction.hh"

JLabKActionInitialization::JLabKActionInitialization()
{}

JLabKActionInitialization::~JLabKActionInitialization()
{}

void JLabKActionInitialization::BuildForMaster() const
{
  SetUserAction(new JLabKRunAction);
}

void JLabKActionInitialization::Build() const
{
  JLabKSteppingAction* steppingAction
  = new JLabKSteppingAction();

  JLabKEventAction* eventAction
  = new JLabKEventAction(steppingAction);

  SetUserAction(new JLabKPrimaryGeneratorAction);
  SetUserAction(new JLabKRunAction);
  SetUserAction(eventAction);
  SetUserAction(steppingAction);
}
