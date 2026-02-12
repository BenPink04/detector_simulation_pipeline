#ifndef JLabKRunAction_hh
#define JLabKRunAction_hh

#include "G4UserRunAction.hh"

class JLabKRunAction : public G4UserRunAction
{
public:
  
  JLabKRunAction();
  virtual ~JLabKRunAction();

  virtual void BeginOfRunAction(const G4Run*);
  virtual void   EndOfRunAction(const G4Run*);

};

#endif
