#include <Rcpp.h>
using namespace Rcpp;



// [[Rcpp::export]]
NumericVector makeInteractionVector(NumericVector x1, NumericVector x2)
{
	int n = x1.size(), i;
	NumericVector x12(n);

	for(i = 0; i < n; i++)
	{
		// x1[i] = x1[i] < 1.0 ? NA_REAL : x1[i] - 1;
		// x2[i] = x2[i] < 1.0 ? NA_REAL : x2[i] - 1;
		x12[i] = x1[i] == 0.0 | x2[i] == 0.0 ? NA_REAL : (x1[i] - 1) * 3 + x2[i] - 1;
		// x12[i] = x1[i] * 3 + x2[i];
	}
	return x12;
}


// [[Rcpp::export]]
NumericMatrix makeCleanYX(NumericVector y, NumericVector x)
{
	int n = x.size(), i, j = 0, nacount = 0;
	NumericVector y1(n), x1(n);

	j = 0;
	for(i = 0; i < n; i++)
	{
		if(NumericVector::is_na(x[i]) | NumericVector::is_na(y[i]))
		{
			nacount++;
		} else {
			y1(j) = y[i];
			x1(j) = x[i];
			j++;
		}
	}

	NumericMatrix out(j-1, 2);
	for(i = 0; i < j-1; i++)
	{
		out(i, 0) = y1[i];
		out(i, 1) = x1[i];
	}
	return(out);
}


// [[Rcpp::export]]
NumericVector ftest8df(NumericVector y, NumericVector X1, NumericVector X2)
{

	int n = X1.size();
	int n1 = X2.size();

	int i, nfac, factor_count[9];
	double MSB, MSW, SSB = 0, SSW = 0, SSD, mY = 0, mY_fac[9], df1, df2;
	NumericVector X, X12, Y, F(3);
	NumericMatrix cleandat;

	X12 = makeInteractionVector(X1, X2);
	cleandat = makeCleanYX(y, X12);

	Y = cleandat( _, 0);
	X = cleandat( _, 1);
	n = X.size();

	for(i = 0; i < 9; i++)
	{
		factor_count[i] = 0;
		mY_fac[i] = 0;
	}
	
	// calculate class means
	for(i = 0; i < n; i++)
	{
		mY += Y[i];
		factor_count[(int)X[i]]++;
		mY_fac[(int)X[i]] += Y[i];
	}

	mY = mY / n;
	nfac = 0;

	// calculate SSB
	for(i = 0; i < 9; i++)
	{
		if (factor_count[i]>0)
		{
			nfac++;
			mY_fac[i] /= factor_count[i];
			SSB += factor_count[i] * pow(mY_fac[i] - mY, 2);
		}
	}
	
	// calculate SSW
	for (i = 0; i < n; i++)
	{
		SSW += pow(Y[i] - mY, 2);
	}


	df1 = nfac - 1;
	df2 = n - df1;
	SSD = SSW - SSB;

	MSB = SSB / df1;
	MSW = SSD / df2;
	F[0] = df1;
	F[1] = df2;
	F[2] = MSB / MSW;

	return F;
}

// [[Rcpp::export]]
NumericVector ftest4df(NumericVector y, NumericVector X1, NumericVector X2)
{
	int n = X1.size();
	int n1 = X2.size();

	int i, nfac, factor_count[9];
	double MSB, MSW, SSB = 0, SSW = 0, SSI = 0, SSD, mY = 0, df1, df2;
	double mY_fac[9];
	double mean_row[3], mean_col[3];

	NumericVector X, X12, Y, F(4);
	NumericMatrix cleandat;

	X12 = makeInteractionVector(X1, X2);
	cleandat = makeCleanYX(y, X12);

	Y = cleandat( _, 0);
	X = cleandat( _, 1);
	n = X.size();


	for(i = 0; i < 9; i++)
	{
		factor_count[i] = 0;
		mY_fac[i] = 0;
	}
	
	// calculate class means
	
	for(i = 0; i < n; i++)
	{
		mY += Y[i];
		factor_count[(int)X[i]]++;
		mY_fac[(int)X[i]] += Y[i];
	}

	mY = mY / n;


	for(i = 0; i < 3; i++)
	{
		mean_col[i] = (mY_fac[i]+mY_fac[i+3]+mY_fac[i+6]) / (factor_count[i]+factor_count[i+3]+factor_count[i+6]);
		mean_row[i] = (mY_fac[3*i]+mY_fac[3*i+1]+mY_fac[3*i+2]) / (factor_count[3*i]+factor_count[3*i+1]+factor_count[3*i+2]);
	}


	// calculate SSB
	nfac = 0;
	for(i = 0; i < 9; i++)
	{
		if(factor_count[i] > 0)
		{
			nfac++;
			mY_fac[i] /= factor_count[i];
			SSB += factor_count[i] * pow(mY_fac[i] - mY, 2);
			SSI += factor_count[i] * pow(mY_fac[i] - mean_row[(int)(i/3)] - mean_col[i%3] + mY, 2);
		}
	}

	// calculate SSW
	for (i = 0; i < n; i++)
	{
		SSW += pow(Y[i] - mY, 2);
	}

	df1 = nfac - 1;
	df2 = n - df1;
	SSD = SSW - SSB;

	MSB = SSB / df1;
	MSW = SSD / df2;

	F[0] = df1;
	F[1] = df2;
	F[2] = MSB / MSW;
	F[3] = (SSI/4) / MSW;

	return F;
}


// [[Rcpp::export]]
NumericMatrix scan8df(NumericVector Y, NumericMatrix X, NumericVector x1)
{
	int p = X.ncol();
	NumericVector x, out;
	NumericMatrix F(p, 3);
	int i;

	for(i = 0; i < p; i++)
	{
		x = X( _, i);
		out = ftest8df(Y, x, x1);
		F(i, 0) = out[0];
		F(i, 1) = out[1];
		F(i, 2) = out[2];
	}

	return F;
}


// [[Rcpp::export]]
NumericMatrix scan4df(NumericVector Y, NumericMatrix X, NumericVector x1)
{
	int p = X.ncol();
	NumericVector x, out;
	NumericMatrix F(p, 4);
	int i;

	for(i = 0; i < p; i++)
	{
		x = X( _, i);
		out = ftest4df(Y, x, x1);
		F(i, 0) = out[0];
		F(i, 1) = out[1];
		F(i, 2) = out[2];
		F(i, 3) = out[3];
	}

	return F;
}

