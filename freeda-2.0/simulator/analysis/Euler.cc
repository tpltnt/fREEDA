#include "Euler.h"

NLEuler::NLEuler(CircVector*& cX, CircVector*& cY, CircVector*& cU,
CircVector*& cTime, double& h)
{
  this->cX = cX;
  this->cY = cY;
  this->cU = cU;
  this->cTime=cTime;
  a = one / h;
  a2 = a * a;
}

NLEuler::NLEuler(CircVector*& cX, CircVector*& cY,
CircVector*& cTime, double& h)
{
  this->cX = cX;
  this->cY = cY;
  this->cTime=cTime;
  a = one / h;
  a2 = a * a;
}

NLEuler::NLEuler(CircVector*& cX, CircVector*& cY, double& h)
{
  //used in SVTran2
  this->cX = cX;
  this->cY = cY;
  this->cTime=NULL;
  a = one / h;
  a2 = a * a;
}

double NLEuler::derivX(const int& index)
{
  return (cX->getCurrent()[index] - cX->getPrevious(1)[index]) * a;
}

double NLEuler::deriv2X(const int& index)
{
  return (cX->getCurrent()[index] - 2. * cX->getPrevious(1)[index]
	+ cX->getPrevious(2)[index]) * a2;
}

double NLEuler::getdx_dtFactor()
{
  return a;
}

double NLEuler::getd2x_dt2Factor()
{
  return a2;
}

double NLEuler::delayX(const int& index, const double& t)
{
  if(cTime==NULL)
	{
		// First, find where to interpolate
		double delta_step = t * a;
		int time_idx = int(delta_step);
		delta_step -= double(time_idx);

		const double& v2 = cX->getPrevious(time_idx)[index];
		const double& v1 = cX->getPrevious(time_idx + 1)[index];
		// Interpolate between v1 and v2
		return v2 + (v1 - v2) * delta_step;
	}
  else
	{
		double curr_time = cTime-> getCurrent()[0];
		int current_pos = cTime->curr_pos;
		double rtime = curr_time-t;
		//define some temporary variables
		double temp ;
		int time_idx;
		int  i=1;
		do
		{
			temp=cTime->getPrevious(i)[0];
			if(curr_time-temp<t)
				// if(temp-time_delay>0)
	    // prev_temp=temp;
	    i++;
		}
		while(curr_time-temp<t && (i<current_pos));
/*
Made fixes here.  i<current_pos above, and time_idx index changes below.
*/

		// while(temp-time_delay>0);
		time_idx = i-1;
		double ttime2 = cTime->getPrevious(time_idx)[0];
		double ttime1 = cTime->getPrevious(time_idx+1)[0];
		double delta_step= abs(rtime-ttime2);
		delta_step =(delta_step)/(ttime2-ttime1);
		// double delta_step=t/(temp-prev_temp);
		//delta_step -=int(delta_step);
		const double& v2 = cX->getPrevious(time_idx)[index];
		const double& v1 = cX->getPrevious(time_idx+1)[index];
		// Interpolate between v1 and v2
		return v2 + (v1 - v2) * delta_step;
	}
}

double NLEuler::getDelayXFactor(const double& t)
{
  if(cTime==NULL)
	{
		// First, find where to interpolate
		double delta_step = t * a;
		int time_idx = int(delta_step);
		delta_step -= double(time_idx);
		if (time_idx)
			return zero;
		else
			return one - delta_step;
	}
  else
	{
		// First, find where to interpolate
		double curr_time = cTime->getCurrent()[0] ;
                int current_pos = cTime->curr_pos;
                double rtime = curr_time-t;
		//define some temporary variables
		double prev_temp, temp ;
		int time_idx;
		int  i=1;
		do
		{
			temp=cTime->getPrevious(i)[0];
			if(curr_time-temp<t)
			  i++;
		}
		while(curr_time-temp<t && (i<current_pos));
		time_idx= i-1;
		double ttime2 = cTime->getPrevious(time_idx)[0];
		double ttime1 = cTime->getPrevious(time_idx+1)[0];
		double delta_step= abs(rtime-ttime2);
		delta_step =(delta_step)/(ttime2-ttime1);
		if (time_idx)
			return zero;
		else
			return one- delta_step;
	}
}

void NLEuler::getNSamples(int &nx, int &ndx)
{
  nx = 1;
  ndx = 0;
}

double NLEuler::derivY(const int& index, const double& currval)
{
  // Save Current value
  cY->getCurrent()[index] = currval;
  return (currval - cY->getPrevious(1)[index]) * a;
}

double NLEuler::deriv(double* x, double * dx)
{
  return a * (x[0] - x[1]);
}

void NLEuler::changeStep(const double& h)
{
  a = one / h;
  a2 = a * a;
}


void NLEuler::predictX(int n_states, double *Integ_predX,
CircVector*& cTimestep)
{
  for(int i=0; i <n_states; i++)
	{
		double temp= (cX->getPrevious(1)[i] - cX->getPrevious(2)[i]);
		double ttime = cTimestep->getPrevious(1)[0];
		temp = temp/ttime;
		//		  /cTimestep->getPrevious(1)[0];
		Integ_predX[i] = cX->getCurrent()[i] = temp * cTimestep->getCurrent()[0] +
		cX->getPrevious(1)[i];
	}
}

void NLEuler::predictU(int ls_size, double *Integ_predU,
CircVector*& cTimestep)
{
  for(int i=0; i <ls_size; i++)
	{
		double temp= (cU->getPrevious(1)[i] - cU->getPrevious(2)[i]);
		temp = temp/cTimestep->getPrevious(1)[0];
		Integ_predU[i] = cU->getCurrent()[i] = temp * cTimestep->getCurrent()[0] +
		cU->getPrevious(1)[i];
	}
}

LEuler::LEuler(CircVector*& cU, TimeMNAM* mnam)
{
  this->cU = cU;
  this->mnam = mnam;
  size = mnam->getDim();
  s2 = DenseDoubleVector(size);
  SerialCommunicator EulerComm;
  ProcessorMap Map(size, 0, EulerComm);
  M1 = new DoubleSparseColMatrix(Copy, Map, int(size/2));
  M1p = new DoubleSparseColMatrix(Copy, Map, int(size/2));
  a = zero;
  rowValExtract = new double[size];
  colIndExtract = new int[size];
}

LEuler::~LEuler()
{
  delete[] rowValExtract;
  delete[] colIndExtract;
  delete M1;
  delete M1p;
}

void LEuler::buildMd(DoubleSparseColMatrix& M, const double& h)
{
  a = one/h;
  mnam->getMatrices(M1, M1p);
  int nnz = 0;
  // copy M1 into M
  for (int i = 0; i < size; i++)
  {
    // extract M1 row-wise
    M1->ExtractGlobalRowCopy(i, size, nnz, rowValExtract, colIndExtract);
    // and put it into M
    M.InsertGlobalValues(i, nnz, rowValExtract, colIndExtract);
  }
  
  // get each row of M1p, scale it by "a" and add it to M
  for (int i = 0; i < size; i++)
  {
    M1p->ExtractGlobalRowCopy(i, size, nnz, rowValExtract, colIndExtract);
    for (int k = 0; k < nnz; k++)
      rowValExtract[k] *= a;
    // finally, put it into M
    // InsertGlobalValues automatically adds to the existing 
    // values in M, but the final FillComplete() call is *required*
    M.InsertGlobalValues(i, nnz, rowValExtract, colIndExtract);
  }
}

void LEuler::buildSf(DenseDoubleVector& s1, const double& ctime)
{
  assert(a);
  mnam->getSource(ctime, s2);
  DistributedDoubleVector u_n(Copy, M1p->RowMap(), &(cU->getPrevious(1)[0]));
  
  // scale u_n by a (wish there was a better way)
  for (int i = 0; i < size; i++)
    u_n[i] *= a; 

  // multiply M1p with u_n and put the result in u_n
  M1p->Multiply(false, u_n, u_n);

  // now add s1 to u_n and put into s1
  for (int i = 0; i < size; ++i)
    s1[i] = s2[i] + u_n[i];
}

void LEuler::buildSf(double* s1, const double& ctime)
{
  assert(a);
  mnam->getSource(ctime, s2);
  memcpy( s1,&(s2[0]),sizeof(double)*size);
  DenseDoubleVector& u_n = cU->getPrevious(1);
  multiply(u_n, s1);
}

void LEuler::multiply(DenseDoubleVector& u_n, double* s1)
{
  col_count=0,row_index=0;
  temp=0;
  Sllist tmp_ptr=NULL, prev_ptr=NULL;
  tmp_ptr = mnam->Lhead;
  prev_ptr=tmp_ptr;

  while(tmp_ptr)
  {
    if(tmp_ptr->row==prev_ptr->row)
    {
      col_count=tmp_ptr->col;
      temp += tmp_ptr->val * u_n[col_count] * a;
      prev_ptr=tmp_ptr;
      tmp_ptr=tmp_ptr->next;
    }
    else
		{
      s1[prev_ptr->row]+=temp;
      prev_ptr=prev_ptr->next;
      temp=0;
    }
  }
  s1[prev_ptr->row]+=temp;
}


void LEuler::changeStep(const double& h)
{
  a = one / h;
}
