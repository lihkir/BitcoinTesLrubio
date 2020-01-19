#include "lminmax_charspeed.h"

void lminmax_charspeed(Matrix<double>& phil, Matrix<double>& phir, Matrix<double>& S)
{
  struct test_cases* pt_test = get_tests();
	
	Matrix<double> &delta = *pt_test->delta; 
	int N = delta.rows();
	
	double phi0 = ((pt_test->nexp - 2)*pt_test->phimax - 1)/(pt_test->nexp - 3);

	double phi1 = sum1D(phil);
	double eta1 = dot_product(delta, phil);
	
	double phi2 = sum1D(phir);
	double eta2 = dot_product(delta, phir);

	double dphi = phi1 - phi2;
	double deta = eta1-eta2;
	
	double w1 = sed_hsf(phi1);
	double w2 = sed_hsf(phi2);

	int i_buf = 1;

	Vector<double> buf_opt(10);

	buf_opt(i_buf) = w1*(delta(1)-eta1); ++i_buf;
	buf_opt(i_buf) = w2*(delta(1)-eta2); ++i_buf;

  double alpha, phi, eta, A, B, C, D, w, min_n, max_n, min_0, max_0;

	if (dphi != 0 && deta != 0)
  {
    /** Hay minimo local en segmento con phi < phi0 ? **/
    alpha = ((pt_test->nexp - 1)*(delta(1) - eta2)*dphi+(1 - phi2)*deta)/(pt_test->nexp*dphi*deta);
    if (alpha > 0 && alpha < 1)
    {
      phi = phi2 + alpha*dphi;
      if (phi < phi0)
      {
        eta = eta2 + alpha*deta;
        buf_opt(i_buf) = sed_hsf(phi)*(delta(1) - eta); ++i_buf;
      }
    }
    
    A = 3*dphi*dphi*deta;
    B = -(2*dphi*deta*(1 + pt_test->phimax - 2*phi2) + 2*dphi*dphi*(delta(1) - eta2));
    C = dphi*(1 + pt_test->phimax - 2*phi2)*(delta(1) - eta2)-(phi2 - pt_test->phimax)*(1 - phi2)*deta;
    /** Hay minimo local en segmento con phi >= phi0 ? **/
    D = B*B-4*A*C;

    if (D >= 0)
    {
      alpha = (-B + sqrt(D))/(2*A);
      if (alpha > 0 && alpha < 1)
      {
        phi = phi2 + alpha*dphi;
        if (phi >= phi0)
        {
          eta = eta2 + alpha*deta;
          buf_opt(i_buf) = sed_hsf(phi)*(delta(1) - eta); ++i_buf;
        }
      }
      
      alpha = (-B - sqrt(D))/(2*A);
      if (alpha > 0 && alpha < 1)
      {
        phi = phi2 + alpha*dphi;
        if (phi >= phi0)
        {
          eta = eta2 + alpha*deta;
          buf_opt(i_buf)=sed_hsf(phi)*(delta(1)-eta);  i_buf=i_buf+1;
        }
      }
    }
  } 
  
  /** singularidad en segmento ? **/
  alpha = -(phi2 - phi0)/dphi;
  if (alpha > 0 && alpha < 1)
  {
    phi = phi2 + alpha*dphi;
    eta = eta2 + alpha*deta;
    buf_opt(i_buf)=sed_hsf(phi)*(delta(1)-eta);  i_buf=i_buf+1;
  }
  
  find_extrema(buf_opt, i_buf - 1, min_0, max_0);

  i_buf = 1;
  /** Extremos segmento **/
  buf_opt(i_buf) = w1*delta(N) + eta1*(sed_hsf_der(phi1)*(1 - phi1) - 2*w1); ++i_buf;
  buf_opt(i_buf) = w2*delta(N) + eta2*(sed_hsf_der(phi2)*(1 - phi2) - 2*w2); ++i_buf;

  if (dphi != 0 && deta != 0)
  {
    /** Hay minimo local en segmento con phi < phi0 ? **/
    alpha = ((pt_test->nexp - 1)*(delta(N) - (pt_test->nexp + 1)*eta2)*dphi + (1 - phi2)*(pt_test->nexp + 1)*deta)/(pt_test->nexp*(pt_test->nexp + 1)*dphi*deta);
    if (alpha > 0 && alpha < 1)
    {
      phi = phi2 + alpha*dphi;
      if (phi < phi0)
      {
        eta = eta2 + alpha*deta;
        w = sed_hsf(phi);
        buf_opt(i_buf)=w*delta(N)+eta*(sed_hsf_der(phi)*(1-phi)-2*w); ++i_buf;
      }
    }

    /** Hay minimo local en segmento con phi >= phi0 ? **/
    C = deta *(1- 5*phi2+ 3*pt_test->phimax + 4*phi2*phi2 - 3*phi2*pt_test->phimax) + dphi*(delta(N)- 5*eta2- 2*delta(N)*phi2+ delta(N)*pt_test->phimax + 8*eta2*phi2 - 3*eta2*pt_test->phimax); 
    B = 2*dphi*(dphi*(4*eta2 - delta(N)) + deta*(-5 + 8*phi2 - 3*pt_test->phimax));
    A = 12*deta*dphi*dphi;
    D = B*B - 4*A*C;
      
    if (D >= 0)
    {
      alpha=(-B + sqrt(D))/(2*A);
      if (alpha > 0 && alpha < 1)
      {
        phi = phi2 + alpha*dphi;
        if (phi >= phi0)
        {
          eta = eta2 + alpha*deta;
          w = sed_hsf(phi);
	        buf_opt(i_buf)=w*delta(N)+eta*(sed_hsf_der(phi)*(1-phi)-2*w); ++i_buf;
        }
      }
        
      alpha = (-B - sqrt(D))/(2*A);
      if (alpha > 0 && alpha < 1)
      {
        phi = phi2 + alpha*dphi;
        if (phi >= phi0)
        {
          eta = eta2 + alpha*deta;
          w = sed_hsf(phi);
          buf_opt(i_buf) = w*delta(N) + eta*(sed_hsf_der(phi)*(1 - phi) - 2*w); ++i_buf;
        }
      }
    }
  }
  /** singularidad en segmento ?**/ 
  alpha = -(phi2 - phi0)/dphi;
  if (alpha > 0 && alpha < 1)
  {
    phi = phi2 + alpha*dphi;
    eta = eta2 + alpha*deta;
    w = sed_hsf(phi);
    buf_opt(i_buf) = w*delta(N) + eta*(sed_hsf_der(phi)*(1 - phi) - 2*w); ++i_buf;
  }
  
  find_extrema(buf_opt, i_buf - 1, min_n, max_n);
  
  S(1) = min_n;
  S(2) = max_0;

  delete pt_test->delta;
	delete pt_test;
}

void find_extrema(Vector<double> &v, int l, double &vmin, double &vmax)
{
  vmin = v(1); vmax = v(1);
  for (int i = 2; i <= l; i++)
  {
    vmin = min(vmin, v(i));
    vmax = max(vmax, v(i));
  }
}