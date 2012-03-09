/***********************************************************************
* State-variable-based Harmonic Balance Analysis Physical Routines
*
*  Author:
*         Carlos E. Christoffersen
***********************************************************************/

#include "SVHB.h"
#include "FreqDomainSV.h"

// Use this for backward compatibility with old routines.
extern "C"
{
	#include "../inout/ftvec.h"
}


/*
* Allocate and fill the T matrix
*/
int SVHB::createT()
{
  int i, j, refcolumn, column;

  // Get total number of states and fill element vector
  int current_ns;
  int n_states = 0;
  max_n_states = 0;
  ElemFlag mask = NONLINEAR;
  // Fill the element vector.
  elem_vec.clear();
  // Reserve a reasonable amount of memory.
  elem_vec.reserve(my_cir->getNumberOfElements()/2);
  // Loop throw all selected elements in circuit.
  my_cir->setFirstElement(mask);
  Element* elem = my_cir->nextElement();
  while(elem)
	{
    // Ask the number of states reported by the element
    n_states += (current_ns = elem->getNumberOfStates());
    // Find the maximum number of states
    if (int(current_ns) > max_n_states)
      max_n_states = current_ns;
    // Put the element in the vector
    elem_vec.push_back(elem);
    // get next element pointer
    elem = my_cir->nextElement();
  }

  // Create interface object
  DenseDoubleVector ov(Teuchos::View, svHBdata->Omega, svHBdata->NoSamples/2);
  fdsv = new FreqDomainSV(ov, svHBdata->NoFreqPoints, max_n_states);

  /* Allocate matrix and fill the matrix with zeros.
	*/
  svHBdata->T = Mlib_CNewMat(n_states, svHBdata->mnam->getDim());
  for (i=0; i < n_states; i++)
	{
    for (j=0; j < int(svHBdata->mnam->getDim()); j++)
		{
      svHBdata->T[i][j].re = svHBdata->T[i][j].im = zero;
		}
	}

  /* Go through all the nonlinear elements and fill the corresponding
	* ones and minus ones.
	*/
  // Number of elements
  unsigned n_elem = elem_vec.size();

  i = 0;
  for (unsigned k=0; k < n_elem; k++)
	{
    UnsignedVector local_ref_vec;
    TerminalVector term_list;
    // Get local reference node information
    elem_vec[k]->getLocalRefIdx(local_ref_vec, term_list);
    // This is just to make sure that things are consistent.
    assert(term_list.size() == elem_vec[k]->getNumTerms());
    // Number of terminal groups (local reference nodes).
    unsigned ngroups = local_ref_vec.size();
    // jbase is the first terminal index in each group
    unsigned jbase = 0;
    for (unsigned l=0; l<ngroups; l++)
		{
      // refcolumn is the local reference terminal column
      refcolumn = term_list[local_ref_vec[l]]->getRC();
      assert(unsigned(refcolumn) <= svHBdata->mnam->getDim());
      for (j = jbase; j < int(local_ref_vec[l]); j++)
			{
				/* Add two matrix elements for each terminal: One element is a one
				* that goes in the terminal column, and the other is a minus one
				* in the local reference terminal column.
				*/
				column = term_list[j]->getRC();
				assert(unsigned(column) <= svHBdata->mnam->getDim());
				assert(i < n_states);
				/* Check that the terminal is not connected
				* to ground before adding the T entry
				*/
				// Remember that MNAM index begin at 0
				if (column)
					svHBdata->T[i][column-1].re = one;
				if (refcolumn)
					svHBdata->T[i][refcolumn-1].re = - one;
				// change state, so we go to the next row.
				i++;
      }
      jbase = local_ref_vec[l] + 1;
    }
  }

  return(n_states);
}

void SVHB::destroyT()
{
  delete fdsv;
  Mlib_CFreeMat(svHBdata->T);
}

void SVHB::get_U_and_I(doublev_t X, dcxm_t *U_NL, dcxm_t *I_NL, dcxm_t *U_L)
{
  int i,j,k;
  int current_ns;

  // Number of elements
  unsigned n_elem = elem_vec.size();

  /* Go through all the nonlinear elements and fill the corresponding
	* U_NL and I_NL matrix elements.
	*/
  i = 0;
  for (unsigned l=0; l < n_elem; l++)
	{
    current_ns = elem_vec[l]->getNumberOfStates();
    for(j = 0; j < current_ns; j++)
		{
      /* put all the frequency components of the j^th state into x[j].
			* x is the time-domain state variable matrix. Each row of the
			* matrix holds a collection of time samples of the values of
			* one state variable. The real parts of the DC component and
			* the last component comes first.
			*/
      fdsv->getX(j)[0] = X[i+j];
      /* Set the value of the last phasor to zero */
      fdsv->getX(j)[1] = zero;
      for (k=1; k < svHBdata->NoFreqPoints; k++)
			{
				fdsv->getX(j)[2*k] = X[(svHBdata->n_states)*(2*k-1) + 2*(i+j)];
				fdsv->getX(j)[2*k+1] = X[(svHBdata->n_states)*(2*k-1) + 2*(i+j)+1];
      }

      /* Now set to zero the rest of the vector. That is
			* because of the oversample we may have more phasors.
			*/
      if (2 * svHBdata->NoFreqPoints < svHBdata->NoSamples)
			{
				for (k= 2 * svHBdata->NoFreqPoints; k < svHBdata->NoSamples; k++)
				{
					fdsv->getX(j)[k] = zero;
				}
			}
    }

    /* Call the svHB evaluation routine of the element.
		*/
    elem_vec[l]->svHB(fdsv);

    /* Fill the frequency-domain vectors.  We save only NoFreqPoints
		* phasors.
		*/
    for(j=0; j< current_ns; j++)
		{
      for (k=0; k < svHBdata->NoFreqPoints; k++)
			{
				svHBws->U_NL[k][i+j].re = fdsv->getVp(j)[2*k];
				svHBws->U_NL[k][i+j].im = fdsv->getVp(j)[2*k+1];

				/* We invert the sign of the currents, in this way the
				* element models can use pasive convention.
				*/
				svHBws->I_NL[k][i+j].re = -fdsv->getIp(j)[2*k];
				svHBws->I_NL[k][i+j].im = -fdsv->getIp(j)[2*k+1];
      }
      /* set the imaginary value of DC to zero */
      svHBws->U_NL[0][i+j].im = zero;
      svHBws->I_NL[0][i+j].im = zero;
    }
    i += current_ns;
  }

  /* Point the U_L, U_NL and I_NL to the matrices in the structure.
	*/
  *U_NL = svHBws->U_NL;
  *I_NL = svHBws->I_NL;
  *U_L = svHBws->U_L;
}

void SVHB::jacobian(double * X, DoubleSparseColMatrix& Jacobian)
{
  DoubleDenseMatrix Jac(svHBdata->sysdim, svHBdata->sysdim);
	// Print an "o" to indicate a Jacobian is being calculated.
  printf("o"); fflush(stdout);

  // Set zero time
  float t0 = (float)clock()/CLOCKS_PER_SEC;

  int nzelem = 0;
  // This two vectors must be cleared
  for (int i=0; i< svHBdata->NoFreqPoints; i++)
    for (int j=0; j< svHBdata->n_states; j++)
      svHBws->JvecU[i][j].re =
			svHBws->JvecU[i][j].im =
			svHBws->JvecI[i][j].re =
			svHBws->JvecI[i][j].im = zero;

	// Number of elements
  unsigned n_elem = elem_vec.size();

  /* Go through all the nonlinear elements and fill the corresponding
	* U_NL and I_NL matrix elements.
	*/
  int i = 0;
  for (unsigned l=0; l < n_elem; l++)
	{
    int current_ns = elem_vec[l]->getNumberOfStates();
    for(int j = 0; j < current_ns; j++)
		{
      /* put all the frequency components of the j^th state into x[j].
			* x is the time-domain state variable matrix. Each row of the
			* matrix holds a collection of time samples of the values of
			* one state variable. The real parts of the DC component and
			* the last component comes first.
			*/
      fdsv->getX(j)[0] = X[i+j];
      /* Set the value of the last phasor to zero */
      fdsv->getX(j)[1] = zero;
      for (int k=1; k < svHBdata->NoFreqPoints; k++)
			{
				fdsv->getX(j)[2*k] = X[(svHBdata->n_states)*(2*k-1) + 2*(i+j)];
				fdsv->getX(j)[2*k+1] = X[(svHBdata->n_states)*(2*k-1) + 2*(i+j)+1];
      }
      /* Now set to zero the rest of the vector. That is
			* because of the oversample we may have more phasors.
			*/
      if (2 * svHBdata->NoFreqPoints < svHBdata->NoSamples)
				for (int k= 2 * svHBdata->NoFreqPoints; k < svHBdata->NoSamples; k++)
					fdsv->getX(j)[k] = zero;
    }

    // Call the deriv_svHB evaluation routine of the element.
    elem_vec[l]->deriv_svHB(fdsv);

    // Fill the Jacobian columns for this element
    for (int jc = 0; jc < current_ns; jc++)
		{
      int j = i + jc;
      for (int findexc = 0; findexc < svHBdata->NoFreqPoints; findexc++)
			{
				for (int m=0; m<2; m++)
				{
	  			if (findexc || !m)
					{
	    			// loop for rows now
	    			for (int jr = 0; jr < current_ns; jr++)
						{
	      			for (int findexr = 0; findexr < svHBdata->NoFreqPoints; findexr++)
							{
								svHBws->JvecU[findexr][i+jr].re =
		  					fdsv->getJacVp(jr, jc)(2*findexr, 2*findexc + m);
								svHBws->JvecU[findexr][i+jr].im =
		  					fdsv->getJacVp(jr, jc)(2*findexr+1, 2*findexc + m);

								svHBws->JvecI[findexr][i+jr].re =
		  					- fdsv->getJacIp(jr, jc)(2*findexr, 2*findexc + m);
								svHBws->JvecI[findexr][i+jr].im =
		  					- fdsv->getJacIp(jr, jc)(2*findexr+1, 2*findexc + m);
	      			}
	    			}
	    			// Calculate Jacobian column here
	    			/* findexr is the freq. index columnwise in the Jacobian.
	     			* i is the state index columnwise.
	     			*/
	    			for (int findexr=0; findexr < svHBdata->NoFreqPoints; findexr++)
						{
	      			/* Now premultiply the current derivative vector by the
	       			* compressed MNA. Instead of using the library function,
	       			* we prefer to use this custom loop that multiply only the
	       			* non-zero elements.
	       			*/
	      			for (int k = 0; k < svHBdata->n_states; k++)
							{
								svHBws->Jtmp[findexr][k].re =
		  					svHBws->Jtmp[findexr][k].im = zero;
								for (int l1 = i; l1 < i+current_ns; l1++)
								{
		  						svHBws->Jtmp[findexr][k].re +=
		    					svHBdata->svMatrix[findexr][k][l1].re *
		    					svHBws->JvecI[findexr][l1].re -
		    					svHBdata->svMatrix[findexr][k][l1].im *
		    					svHBws->JvecI[findexr][l1].im;
		  						svHBws->Jtmp[findexr][k].im +=
		    					svHBdata->svMatrix[findexr][k][l1].re *
		    					svHBws->JvecI[findexr][l1].im +
		    					svHBdata->svMatrix[findexr][k][l1].im *
		    					svHBws->JvecI[findexr][l1].re;
								}
	      			}

	      			// Now fill columns of the Jacobian.
	      			for (int i1=0; i1< svHBdata->n_states; i1++)
							{
								int column = (!findexc) ?
		  					j : (2*findexc-1) * (svHBdata->n_states) + 2*j + m;
								int row = (!findexr) ?
		  					i1 : (2*findexr-1) * (svHBdata->n_states) + 2*i1;

								if (fabs(Jac(row,column) =
									svHBws->Jtmp[findexr][i1].re
			 					- svHBws->JvecU[findexr][i1].re) >
		    				1.e-3)
								nzelem++;

								if (findexr)
		  						if (fabs(Jac(row+1, column) =
			   						svHBws->Jtmp[findexr][i1].im
								- svHBws->JvecU[findexr][i1].im) >
								1.e-3)
								nzelem++;
	      			}
	    			}
	  			}
				}
      }
    }
    // Clean vectors
    for (int jr = 0; jr < current_ns; jr++)
		{
      for (int findexr = 0; findexr < svHBdata->NoFreqPoints; findexr++)
			{
				svHBws->JvecU[findexr][i+jr].re =
				svHBws->JvecU[findexr][i+jr].im =
				svHBws->JvecI[findexr][i+jr].re =
				svHBws->JvecI[findexr][i+jr].im = zero;
      }
    }
    // Increment the state variable counter
    i += current_ns;
  }

  float density = (float) nzelem / (svHBdata->sysdim * svHBdata->sysdim);

  // copy contents of dense jacobian J into a sparse equivalent Jacobian
  // required by NOX
  std::vector<int> indices(svHBdata->sysdim, 0);
  std::vector<double> values(svHBdata->sysdim, 0.0);
  
  int k;
  // row loop
  for (int i = 0; i < svHBdata->sysdim; i++)
  {
    k = 0;
    // column loop
    for(int j = 0; j < svHBdata->sysdim; j++)
    {
      if (Jac(i,j) != 0)
      {
        // We have a non-zero value save it
        values[k] = Jac(i,j);
        indices[k] = j; //Save the column the value was found in.
        k++;
      }
    }
    if (k != 0)
    {
      // Non zero values were found for row i
      Jacobian.ReplaceMyValues(i, k, &values[0], &indices[0]);
    }
  }
  if (svHBopt->verbosity > 0)
  {
    fprintf(output_F,"   * Jacobian Density: %3.0f %%, CPU Time: %6.2f s\n",
        density*100, (float)clock()/CLOCKS_PER_SEC - t0);
    fflush(output_F);
  }
}

void SVHB::doOutput(dcxm_t VI, dcxm_t I_NL,
doublev_t X, double residual)
{
  unsigned i, j;
  int findex;
  int NoFreq = 0;
  unsigned order;

  /* If NoFreqPoints is different than 2^n+1 and we did a single-tone
	* analysis, then make some special processing to keep the output
	* invfft happy.
	*/
  if (svHBopt->numtones == 1)
	{
    NoFreq = 0;
    /* Search the minimum power of 2 greater than noFreqPoints */
    for (order=1; NoFreq < svHBopt->h[0]+1; order++)
      NoFreq = (1 << order);
    NoFreq++;
    if(svHBdata->NoFreqPoints < NoFreq)
		{
      /* allocate memory for the frequency vector
			*/
      allocFreqV_P(NoFreq);
      /* fill vector */
      for (findex=0; findex < NoFreq; findex++)
				FreqV_P[findex] = svHBopt->tone[0] * findex;
    }
  }
  else
	{
    NoFreq = svHBdata->NoFreqPoints;
    /* allocate memory for the frequency vector
		*/
    allocFreqV_P(svHBdata->NoFreqPoints);
    /* fill vector */
    for (findex=0; findex < svHBdata->NoFreqPoints; findex++)
      FreqV_P[findex] = svHBdata->FreqV_P[findex];
  }

  // Temporary vector for voltages (initialized to zero)
  DenseComplexVector tmp_v(NoFreq);

  // For each terminal, assign voltage vector
  Terminal* term = NULL;
  my_cir->setFirstTerminal();
  while((term = my_cir->nextTerminal()))
	{
    // Get MNAM index
    if (term->getRC())
		{
      // Remember that MNAM indices begin at 0
      i = term->getRC() - 1;
      // Get vector from VI
      for (findex=0; findex < svHBdata->NoFreqPoints; findex++)
      {
        tmp_v[findex].real() = VI[findex][i].re;
        tmp_v[findex].imag() = VI[findex][i].im;
      }
    }
    else
      tmp_v.putScalar(double_complex(zero)); // This is a reference terminal
    
    // Set terminal vector
    term->getTermData()->setPhasorV(tmp_v);
  }

  // Temporary vectors (initialized to zero)
  DenseComplexVector tmp_i(NoFreq);
  DenseComplexVector tmp_x(NoFreq);

  ElemFlag mask = LINEAR;
  // Loop throw all selected elements in circuit and fill current
  // vector if needed.
  my_cir->setFirstElement(mask);
  Element* elem = my_cir->nextElement();
  unsigned first_eqn, no_eqn;
  while(elem)
	{
    elem->getExtraRC(first_eqn, no_eqn);
    if (first_eqn)
		{
      // Fill current vector(s) in element
      for (i = 0; i < no_eqn; i++)
			{
				// Remember that MNAM indices begin at 0
				unsigned row = first_eqn - 1 + i;
				for (findex=0; findex < svHBdata->NoFreqPoints; findex++)
        {
          tmp_i[findex].real() = VI[findex][row].re;
          tmp_i[findex].imag() = VI[findex][row].im;
        }
				elem->getElemData()->setPhasorI(i, tmp_i);
      }
    }
    // get next element pointer
    elem = my_cir->nextElement();
  }

  /* Now fill node vectors of nonlinear devices.
	*/
  unsigned current_ns;
  unsigned n_elem = elem_vec.size();   // Number of elements
  i=0;
  for (unsigned k=0; k < n_elem; k++)
	{
    current_ns = elem_vec[k]->getNumberOfStates();
    for(j=0; j< current_ns; j++)
		{
      /* Fill first the DC component */
      tmp_x[0].real() = X[i+j];
      tmp_x[0].imag() = zero;

      /* For the current, invert the sign to keep passive convention.
			*/
      tmp_i[0].real() = -I_NL[0][j+i].re;
      tmp_i[0].imag() = zero;

      /* Now fill the rest of the harmonics */
      for (findex=1; findex < svHBdata->NoFreqPoints; findex++)
			{
        tmp_x[findex].real() = X[(svHBdata->n_states)*(2*findex-1) + 2*(i+j)]; 
        tmp_x[findex].imag() = X[(svHBdata->n_states)*(2*findex-1) + 2*(i+j)+1];

        tmp_i[findex].real() = -I_NL[findex][i+j].re;
        tmp_i[findex].imag() = -I_NL[findex][i+j].im;
      }
      elem_vec[k]->getElemData()->setPhasorI(j, tmp_i);
      elem_vec[k]->getElemData()->setPhasorX(j, tmp_x);
    }
    i += current_ns;
  }

  /* Print the residual */
  printf("\n\n *** Final norm of the residuals: %g\n\n", residual);
  fprintf(output_F, "\n *** Final norm of the residuals: %g\n\n", residual);
}

void SVHB::createWs()
{
  svHBws = (svHB_work_space_t *) malloc(sizeof(svHB_work_space_t));
  /* Allocate memory for frequency-domain vectors.
	*/
  svHBws->JvecU = Mlib_CNewMat(svHBdata->NoFreqPoints, svHBdata->n_states);
  svHBws->JvecI = Mlib_CNewMat(svHBdata->NoFreqPoints, svHBdata->n_states);
  svHBws->Jtmp = Mlib_CNewMat(svHBdata->NoFreqPoints, svHBdata->n_states);

  /* Share the space with the Jacobian vectors, we don't use them
	* at the same time.
	*/
  svHBws->U_NL = svHBws->JvecU;
  svHBws->I_NL = svHBws->JvecI;
  svHBws->U_L = svHBws->Jtmp;
}

void SVHB::destroyWs()
{
  Mlib_CFreeMat(svHBws->JvecU);
  Mlib_CFreeMat(svHBws->JvecI);
  Mlib_CFreeMat(svHBws->Jtmp);
  free(svHBws);
}

