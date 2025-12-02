#ifndef JLabKEventAction_hh
#define JLabKEventAction_hh

#include "G4UserEventAction.hh"
#include "globals.hh"

class JLabKRunAction;
class JLabKVSteppingAction;

class JLabKEventAction : public G4UserEventAction
{
public:

  JLabKEventAction(JLabKVSteppingAction*);

  virtual ~JLabKEventAction();

  virtual void BeginOfEventAction(const G4Event*);
  virtual void EndOfEventAction(const G4Event*);

private:

  JLabKVSteppingAction* fpJLabKVSteppingAction;
};

#endif
