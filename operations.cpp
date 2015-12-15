#include <Rcpp.h>

using namespace Rcpp;

// trace of a matrix
// [[Rcpp::export]]
double trace(NumericMatrix x)
{
    double t = 0;
    int p = x.nrow();
    for (int i=0; i<p; i++)
    {
        t += x(i,i);
    }
    return t;
}


// inner product of two vectors
// [[Rcpp::export]]
double inner(NumericVector x, NumericVector y) 
{
    int nx = x.size();
    int ny = y.size();
    double sumup=0;

    if(nx==ny)
    {
      for(int i=0; i < nx; i++)
          sumup += x[i] * y[i];
    }
    else
        sumup = NA_REAL; // NA_REAL: constant of NA value for numeric (double) values
  
    return sumup;
}


// colMeans in C++
// [[Rcpp::export]]
NumericVector colMeansC(NumericMatrix x) 
{
  int nrow = x.nrow(), ncol = x.ncol();
  NumericVector out(ncol);
  for (int j = 0; j < ncol; j++) 
  {
    double total = 0;
    for (int i = 0; i < nrow; i++) 
    {
      total += x(i, j);
    }
    out[j] = total;
  }
  return out / nrow;
}

// Akr
// [[Rcpp::export]]
double App(NumericVector x, NumericVector y)
{
    NumericVector d = x-y;
    double a = inner(d,d);
    return a;
}

// A2krls
// [[Rcpp::export]]
double A2pp(NumericVector x, NumericVector y, NumericVector u, NumericVector v)
{
    NumericVector d1 = x-y;
    NumericVector d2 = u-v;
    double a = pow(inner(d1,d2),2);
    return a;
}

// B2kr
// [[Rcpp::export]]
double Bpp(NumericVector x, NumericVector y)
{
    double b = pow(inner(x,y),2); 
    return b;
}


// remove rows with indices rowID from a matrix x
// [[Rcpp::export]]
NumericMatrix row_erase(NumericMatrix x, IntegerVector rowID) 
{
  rowID = rowID - 1;
  rowID = rowID.sort();
  NumericMatrix x2(Dimension(x.nrow()- rowID.size(), x.ncol()));
  int iter = 0; 
  int del = 1; // to count deleted elements
  for (int i = 0; i < x.nrow(); i++) {
    if (i != rowID[del - 1]) {
      x2.row(iter) = x.row(i);
      iter++;
    } else {
      del++;
    }
  }
  return x2;
}

// D(k,r)
// [[Rcpp::export]]
NumericVector Dppp(NumericMatrix x, int k, int r)
{
    IntegerVector kr = IntegerVector::create(k,r);
    NumericMatrix y = row_erase(x, kr);
    NumericVector d = x.row(k-1) - colMeansC(y);
    return d;
}


// C2kr
// [[Rcpp::export]]
double Cpp(NumericMatrix x, int k, int r)
{
    double c = inner(x.row(k-1),Dppp(x,k,r))*inner(x.row(r-1),Dppp(x,r,k));
    return c;
}


// B2krls
// [[Rcpp::export]]
double B2pp(NumericVector x, NumericVector y, NumericVector u, NumericVector v)
{
    double b = A2pp(x,y,u,v) + A2pp(x,u,y,v) + A2pp(x,v,u,y);
    return b;
}

// Q(n)
// [[Rcpp::export]]
int Qpp(int ni)
{
    int q = ni*(ni-1);
    return q;
}


// P(n)
// [[Rcpp::export]]
int Ppp(int ni)
{
    int p = Qpp(ni)*(ni-2)*(ni-3);
    return p;
}


// E3
// [[Rcpp::export]]
double E3pp(NumericMatrix x, IntegerMatrix perm) 
{
    int n = x.nrow();
    int nperm = perm.nrow();
    double e3 = 0;
    for (int i = 0; i < nperm; i++)
    {
        int k = perm(i,0);
        int r = perm(i,1);
        NumericVector xk = x.row(k-1);
        NumericVector xr = x.row(r-1);
        e3 += Bpp(xk,xr);
    }
    return e3 / Qpp(n);
}


// E4
// [[Rcpp::export]]
double E4pp(NumericMatrix x, IntegerMatrix perm) 
{
    int n = x.nrow();
    int nperm = perm.nrow();
    double e4 = 0;
    for (int i = 0; i < nperm; i++)
    {
        int k = perm(i,0);
        int r = perm(i,1);
        e4 += Cpp(x,k,r);
    }
    return e4 / Qpp(n);
}

// E5
// [[Rcpp::export]]
double E5pp(NumericMatrix x, IntegerMatrix perm) 
{
    int n = x.nrow();
    int nperm = perm.nrow();
    double e5 = 0;
    for (int i = 0; i < nperm; i++)
    {
        int k = perm(i,0);
        int r = perm(i,1);
        int l = perm(i,2);
        int s = perm(i,3);
        e5 += A2pp(x.row(k-1), x.row(r-1), x.row(l-1), x.row(s-1));
    }
    return e5 / (4*Ppp(n));
}


// E6
// [[Rcpp::export]]
double E6pp(NumericMatrix x, IntegerMatrix perm) 
{
    int n = x.nrow();
    int nperm = perm.nrow();
    double e6 = 0;
    for (int i = 0; i < nperm; i++)
    {
        int k = perm(i,0);
        int r = perm(i,1);
        int l = perm(i,2);
        int s = perm(i,3);
        e6 += B2pp(x.row(k-1), x.row(r-1), x.row(l-1), x.row(s-1));
    }
    return e6 / (12*Ppp(n));
}
