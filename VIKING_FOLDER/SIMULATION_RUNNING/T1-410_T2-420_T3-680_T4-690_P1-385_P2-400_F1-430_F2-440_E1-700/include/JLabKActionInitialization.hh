#ifndef JLabKActionInitialization_hh
#define JLabKActionInitialization_hh 1

#include "G4VUserActionInitialization.hh"

class JLabKActionInitialization : public G4VUserActionInitialization
{
public:
  JLabKActionInitialization();
  virtual ~JLabKActionInitialization();

  virtual void BuildForMaster() const;
  virtual void Build() const;
};

#endif


