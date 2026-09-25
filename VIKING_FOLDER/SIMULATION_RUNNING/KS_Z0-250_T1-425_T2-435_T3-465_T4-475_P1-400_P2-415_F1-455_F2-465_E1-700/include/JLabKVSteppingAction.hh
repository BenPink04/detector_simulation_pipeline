#ifndef JLabKVSteppingAction_hh
#define JLabKVSteppingAction_hh

#include "G4UserSteppingAction.hh"
#include "globals.hh"

class JLabKVSteppingAction : public G4UserSteppingAction
{
public:
  virtual void BeginOfEventAction() {};
  virtual void EndOfEventAction() {};
};

#endif
