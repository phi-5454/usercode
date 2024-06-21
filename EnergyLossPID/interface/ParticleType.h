#ifndef _ParticleType_h_
#define _ParticleType_h_

class MostProbable;

class ParticleType
{
 public:
  ParticleType(MostProbable * mostProbable_, int pass_);
  virtual ~ParticleType();

  int guess(double p, double y);
  int sure (double p, double y);

 private:
  MostProbable * mostProbable;
  int pass;
};

#endif
