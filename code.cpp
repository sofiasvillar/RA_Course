
#include <vector>
#include <Rcpp.h>

// [[Rcpp::plugins("cpp11")]]

enum methodT {
    ER=0,
    PW=1,
    RPW=2,
    BRAR=3,
    SMLE_Neyman=4,
    SMLE_minF=5,
    SMLE_AD=6,
    DBCD_Neyman=7,
    DBCD_minF=8,
    DBCD_AD=9,
    ERADE_Neyman=10,
    ERADE_minF=11,
    ERADE_AD=12
};
    
enum measureT {
    SMD=0,
    LRR=1
};

enum arT {
    NEYMAN=0,
    MINF=1,
    AD=2
};

using std::vector;

const double INF = std::numeric_limits<double>::infinity();
const double NEGINF = -1 * std::numeric_limits<double>::infinity();

// Function superior(p1, p2) = 1 if p1 < p2; 0 if p1 > p2; 2 if equal

double setinf(double p0, double p1) {
    if (p0 == p1)
	return 0;
    else if (p0 < p1)
	return NEGINF;
    else // p0 > p1)
	return INF;
}

double pnorm(double in) {
    return R::pnorm(in, 0.0, 1.0, true, false);
}

double wald(std::vector<double> &x0, std::vector<double> &x1, int measure) {   
   
    // Vector mean
    double p0 = std::accumulate(x0.begin(), x0.end(), 0) / (double) x0.size();
    double p1 = std::accumulate(x1.begin(), x1.end(), 0) / (double) x1.size();
    double q0 = 1 - p0;
    double q1 = 1 - p1;
    double sdx0 = p0 * q0;
    double sdx1 = p1 * q1;
    int n0 = x0.size();
    int n1 = x1.size();
    
    double Z;
    
    if (measure == SMD) {
	if (sdx0 == 0 && sdx1 == 0)
	    Z = setinf(p0, p1);
	else
	    Z = (p0 - p1) / sqrt(sdx0 / n0 + sdx1 / n1);
	    
    } else if (measure == LRR) {
	if (p0 == 1 || p1 == 1 || (p0 == 0 && p1 == 0))
	    Z = setinf(p0, p1);
	else
	    Z = log(q1 / q0) / sqrt(p0/(n0*q0) + p1/(n1*q1));	
	
    } else {
	Rcpp::Rcout << "METHOD NOT SUPPORTED\n";
    }
    return 2 * pnorm( -1 * fabs(Z));
}

double wald(Rcpp::NumericVector &x0, Rcpp::NumericVector &x1, int measure) {
    vector<double> v0 = Rcpp::as<vector<double>>(x0);
    vector<double> v1 = Rcpp::as<vector<double>>(x1);
    return wald(v0, v1, measure);
}

std::vector<double> two_arm_ER(int N, double p1, double p2, int measure, int burnin, bool Zcorr) {

    // rbinom returns double
    int n1 = std::min(std::max((int) std::lround(R::rbinom(N, 0.5)), burnin), N-burnin);
    int n0 = N - n1;
    Rcpp::NumericVector x0 = Rcpp::rbinom(n0, 1, p1);
    Rcpp::NumericVector x1 = Rcpp::rbinom(n1, 1, p2);

    int s0 = 0;
    for (int val : x0)
	s0 += val;
    int s1 = 0;
    for (int val : x1)
	s1 += val;
    double response = (s0 + s1) / (double) N;

    double per_sup;
    if (p1 > p2)
	per_sup = n0 / (double) N;
    else if (p1 < p2)
	per_sup = n1 / (double) N;
    else
	per_sup = NA_REAL;
    
    double Z_P = wald(x0, x1, measure);
    
    vector<double> out = {response, Z_P, per_sup};
    return out;
}


std::vector<double> PWfun(int N, double p1, double p2, int burnin, int measure, bool Zcorr) {
    std::vector<int> allocation(N);
    std::vector<double> response(N);

    // burnin > 0 is enforced
    for (int i = 0; i < burnin; i++)
	response[i] = R::rbinom(1, p1);
    for (int i = burnin; i < 2*burnin; i++) {
	allocation[i] = 1;
	response[i] = R::rbinom(1, p2);
    }
    for (int i = 2*burnin; i < N; i++) {
	if (response[i-1] == 1)
	    allocation[i] = allocation[i-1];
	else
	    allocation[i] = abs(allocation[i-1] - 1);
	if (allocation[i] + 1 == 1)
	    response[i] = R::rbinom(1, p1);
	else if (allocation[i] + 1 == 2)
	    response[i] = R::rbinom(1, p2);
	else
	    Rcpp::Rcout << "ERROR PW rbinom probs\n";
    }

    int n0 = 0, n1 = 0, s0 = 0, s1 = 0;
    std::vector<double> x0, x1;
    for (int i = 0; i < allocation.size(); i++) {
	if (allocation[i] == 0) {
	    n0++;
	    s0 += response[i];
	    x0.push_back(response[i]);
	} else { // allocation == 1
	    n1++;
	    s1 += response[i];
	    x1.push_back(response[i]);
	}
    }
    double response2 = (s0 + s1) / (double) N;
   
    // Function superior(p1, p2) = 1 if p1 < p2; 0 if p1 > p2; 2 if equal

    double per_sup;
    if (p1 > p2)
	per_sup = n0 / (double) N;
    else if (p1 < p2)
	per_sup = n1 / (double) N;
    else
	per_sup = NA_REAL;

    double Z_P = wald(x0, x1, measure);

    vector<double> out = {response2, Z_P, per_sup};
    return out;
    
}


std::vector<double> randomiseRPW(int N, double p1, double p2, int burnin, int measure, bool zcorr) {

    std::vector<int> allocation(N);
    std::vector<double> response(N);
    std::vector<int> urn(N);
    
    for (int i = 0; i < burnin; i++) {
	response[i] = R::rbinom(1, p1);
	urn[i] = 0;
    }
    for (int i = burnin; i < 2*burnin; i++) {
	allocation[i] = 1;
	response[i] = R::rbinom(1, p2);
	urn[i] = 1;
    }
    for (int i = 2*burnin; i < N; i++) {
	// Call R sample to preserve RNG
	auto start = urn.begin();
	auto end = urn.begin() + i;
	Rcpp::IntegerVector subvec(start, end);
	int ball = Rcpp::sample(subvec, 1)[0];
	allocation[i] = ball;
	if (allocation[i] + 1 == 1)
	    response[i] = R::rbinom(1, p1);
	else if (allocation[i] + 1 == 2)
	    response[i] = R::rbinom(1, p2);
	urn[i] = response[i] ? ball : 1-ball;
    }

    int n0 = 0, n1 = 0, s0 = 0, s1 = 0;
    std::vector<double> x0, x1;
    for (int i = 0; i < allocation.size(); i++) {
	if (allocation[i] == 0) {
	    n0++;
	    s0 += response[i];
	    x0.push_back(response[i]);
	} else { // allocation == 1
	    n1++;
	    s1 += response[i];
	    x1.push_back(response[i]);
	}
    }
    double response2 = (s0 + s1) / (double) N;
    
    double per_sup;
    if (p1 > p2)
	per_sup = n0 / (double) N;
    else if (p1 < p2)
	per_sup = n1 / (double) N;
    else
	per_sup = NA_REAL;
    
    double Z_P = wald(x0, x1, measure);

    vector<double> out = {response2, Z_P, per_sup};
    return out;
    
}


int selectArm(int s1, int s2, int f1, int f2) {
    double val1 = R::rbeta(s1+1, f1+1);
    double val2 = R::rbeta(s2+1, f2+1);
    if (val1 >= val2)
	return 1;
    else
	return 2;
}


std::vector<double> BRARfun(int N, double p1, double p2, int burnin, int measure, bool zcorr) {
    vector<int> A(N), X(N);

    for (int i = 0; i < burnin; i++) {
	A[i] = 1;
	X[i] = R::rbinom(1, p1);
    }
    for (int i = burnin; i < 2*burnin; i++) {
	A[i] = 2;
	X[i] = R::rbinom(1, p2);
    }

    for (int i = 2*burnin; i < N; i++) {
	int s1 = 0, s2 = 0, tot1 = 0, tot2 = 0, f1 = 0, f2 = 0;
	for (int j = 0; j < N; j++) {
	    if (A[j] == 1) {
		s1 += X[j];
		tot1++;
	    }
	    else if (A[j] == 2) {
		s2 += X[j];
		tot2++;
	    } else { // (A[j] == 0, equiv to NA and counted!
		tot1++;
		tot2++;
	    }
	}
	f1 = tot1 - s1;
	f2 = tot2 - s2;
	A[i] = selectArm(s1, s2, f1, f2);
	if (A[i] == 1)	
	    X[i] = R::rbinom(1, p1);
	else // ==2
	    X[i] = R::rbinom(1, p2);
    }

    int n0 = 0, n1 = 0, s0 = 0, s1 = 0;
    std::vector<double> x0, x1;
    for (int i = 0; i < A.size(); i++) {
	if (A[i] == 1) {
	    n0++;
	    s0 += X[i];
	    x0.push_back(X[i]);
	} else { // == 2
	    n1++;
	    s1 += X[i];
	    x1.push_back(X[i]);
	}
    }
    double response = (s0 + s1) / (double) N;
    
    double per_sup;
    if (p1 > p2)
	per_sup = n0 / (double) N;
    else if (p1 < p2)
	per_sup = n1 / (double) N;
    else
	per_sup = NA_REAL;
    
    double Z_P = wald(x0, x1, measure);

    vector<double> out = {response, Z_P, per_sup};
    return out;
    
}



vector<double> two_arm_SMLE(int N, double p1, double p2, int measure, int burnin, int ar) {
    vector<int> A(N), X(N);

    for (int i = 0; i < burnin; i++) {
	A[i] = 1;
	X[i] = R::rbinom(1, p1);
    }
    for (int i = burnin; i < 2*burnin; i++) {
	A[i] = 2;
	X[i] = R::rbinom(1, p2);
    }
    
    for (int i = 2*burnin; i < N; i++) {
	int n1 = 0, n2 = 0;
	for (int j = 0; j < A.size(); j++) {
	    if (A[j] == 1)
		n1++;
	    if (A[j] == 2)
		n2++;
	}
	double p;
	if (ar == NEYMAN) {
	    if (measure == SMD) {
		double sig1_N = 0, sig2_N = 0;
		double mean1 = 0, mean2 = 0;
		int tot1 = 0, tot2 = 0;
		for (int j = 0; j < A.size(); j++) {
		    if (A[j] == 1) {
			mean1 += X[j]; tot1++;
		    } else if (A[j] == 2) {
			mean2 += X[j]; tot2++;
		    }
		}
		mean1 /= tot1; mean2 /= tot2;
		for (int j = 0; j < A.size(); j++) {
		    if (A[j] == 1)
			sig1_N += pow(X[j] - mean1, 2);
		    else if (A[j] == 2)
			sig2_N += pow(X[j] - mean2, 2);
		}
		sig1_N = sqrt(sig1_N / (tot1-1));
		sig2_N= sqrt(sig2_N / (tot2-1));
		
		
		if (sig1_N + sig2_N == 0)
		    p = 0.5;
		else
		    p = 1 - sig1_N / (sig1_N + sig2_N);
	    } else {
		int s1 = 0, s2 = 0;
		for (int j = 0; j < A.size(); j++) {
		    if (A[j] == 1)
			s1 += X[j];
		    else if (A[j] == 2)
			s2 += X[j];
		}
		double p1hat = s1 / (double) n1;
		double p2hat = s2 / (double) n2;
		if ( (sqrt((1-p1hat)*p2hat) + sqrt((1-p2hat)*p1hat)) == 0)
		    p = 0.5;
		else
		    p = sqrt((1-p1hat)*p2hat) / (sqrt((1-p1hat)*p2hat) + sqrt((1-p2hat)*p1hat));
	    }
	} else if (ar == MINF) {
	    int s1 = 0, s2 = 0;
	    for (int j = 0; j < A.size(); j++) {
		if (A[j] == 1)
		    s1 += X[j];
		else if (A[j] == 2)
		    s2 += X[j];
	    }
	    double p1hat = s1 / (double) n1;
	    double p2hat = s2 / (double) n2;

	    double rho1hat;
	    if (measure == SMD) {
		if (p1hat + p2hat == 0)
		    rho1hat = 0.5;
		else
		    rho1hat = sqrt(p1hat) / (sqrt(p1hat) + sqrt(p2hat));
	    } else { // LRR
		if ((sqrt(p1hat)*(1-p2hat) + sqrt(p2hat)*(1-p1hat)) == 0)
		    rho1hat = 0.5;
		else
		    rho1hat = 1 - sqrt(p2hat)*(1-p1hat) / (sqrt(p1hat)*(1-p2hat) + sqrt(p2hat)*(1-p1hat));
	    }
	    p = 1 - rho1hat;
	} else { // (ar == AD)
	    if (measure == SMD) {
		int s1 = 0, s2 = 0;
		for (int j = 0; j < A.size(); j++) {
		    if (A[j] == 1)
			s1 += X[j];
		    else if (A[j] == 2)
			s2 += X[j];
		}
		double p1hat = s1 / (double) n1;
		double p2hat = s2 / (double) n2;

		double rho1hat = 0;
		if (p1hat + p2hat == 0)
		    rho1hat = 0.5;
		else
		    rho1hat = p1hat / (p1hat + p2hat);
		p = 1 - rho1hat;

	    } else
		// NOT SUPPORTED
		p = 0;
	}
	
	if (p == 0)
	    p = 1/(double)N;
	if (p == 1)
	    p = 1-1/(double)N;
	A[i] = R::rbinom(1,p)+1;
	if (A[i] == 1)
	    X[i] = R::rbinom(1, p1);
	else
	    X[i] = R::rbinom(1, p2);
    }
    
    int n0 = 0, n1 = 0, s0 = 0, s1 = 0;
    std::vector<double> x0, x1;
    for (int i = 0; i < A.size(); i++) {
	if (A[i] == 1) {
	    n0++;
	    s0 += X[i];
	    x0.push_back(X[i]);
	} else { // == 2
	    n1++;
	    s1 += X[i];
	    x1.push_back(X[i]);
	}
    }
    double response = (s0 + s1) / (double) N;
    
    double per_sup;
    if (p1 > p2)
	per_sup = n0 / (double) N;
    else if (p1 < p2)
	per_sup = n1 / (double) N;
    else
	per_sup = NA_REAL;
    
    double Z_P = wald(x0, x1, measure);
    
    vector<double> out = {response, Z_P, per_sup};
    return out;
 	    
}


vector<double> two_arm_DBCD(int N, double p1, double p2, int measure, int burnin, int ar) {
    vector<int> A(N), X(N);

    for (int i = 0; i < burnin; i++) {
	A[i] = 1;
	X[i] = R::rbinom(1, p1);
    }
    for (int i = burnin; i < 2*burnin; i++) {
	A[i] = 2;
	X[i] = R::rbinom(1, p2);
    }
    
   for (int i = 2*burnin; i < N; i++) {
	int n1 = 0, n2 = 0;
	for (int j = 0; j < A.size(); j++) {
	    if (A[j] == 1)
		n1++;
	    if (A[j] == 2)
		n2++;
	}

	double p;
	
	if (ar == NEYMAN) {
	    if (measure == SMD) {
		double sig1_N = 0, sig2_N = 0;
		double mean1 = 0, mean2 = 0;
		int tot1 = 0, tot2 = 0;
		for (int j = 0; j < A.size(); j++) {
		    if (A[j] == 1) {
			mean1 += X[j]; tot1++;
		    } else if (A[j] == 2) {
			mean2 += X[j]; tot2++;
		    }
		}
		mean1 /= tot1; mean2 /= tot2;
		for (int j = 0; j < A.size(); j++) {
		    if (A[j] == 1)
			sig1_N += pow(X[j] - mean1, 2);
		    else if (A[j] == 2)
			sig2_N += pow(X[j] - mean2, 2);
		}
		sig1_N = sqrt(sig1_N / (tot1-1));
		sig2_N = sqrt(sig2_N / (tot2-1));
		
		if (sig1_N + sig2_N == 0)
		    p = 0.5;
		else
		    p = 1 - sig1_N / (sig1_N + sig2_N);
	    } else { // LRR
		int s1 = 0, s2 = 0;
		for (int j = 0; j < A.size(); j++) {
		    if (A[j] == 1)
			s1 += X[j];
		    else if (A[j] == 2)
			s2 += X[j];
		}
		double p1hat = s1 / (double) n1;
		double p2hat = s2 / (double) n2;
		if ( (sqrt((1-p1hat)*p2hat) + sqrt((1-p2hat)*p1hat)) == 0)
		    p = 0.5;
		else
		    p = sqrt((1-p1hat)*p2hat) / (sqrt((1-p1hat)*p2hat) + sqrt((1-p2hat)*p1hat));
	    }

	} else if (ar == MINF) {
	    int s1 = 0, s2 = 0;
	    for (int j = 0; j < A.size(); j++) {
		if (A[j] == 1)
		    s1 += X[j];
		else if (A[j] == 2)
		    s2 += X[j];
	    }
	    double p1hat = s1 / (double) n1;
	    double p2hat = s2 / (double) n2;

	    double rho1hat;
	    if (measure == SMD) {
		if (p1hat + p2hat == 0)
		    rho1hat = 0.5;
		else
		    rho1hat = sqrt(p1hat) / (sqrt(p1hat) + sqrt(p2hat));
	    } else { // LRR
		if ((sqrt(p1hat)*(1-p2hat) + sqrt(p2hat)*(1-p1hat)) == 0)
		    rho1hat = 0.5;
		else
		    rho1hat = 1 - sqrt(p2hat)*(1-p1hat) / (sqrt(p1hat)*(1-p2hat) + sqrt(p2hat)*(1-p1hat));
	    }
	    p = 1 - rho1hat;

	} else { // (ar == AD)
	    if (measure == SMD) {
		int s1 = 0, s2 = 0;
		for (int j = 0; j < A.size(); j++) {
		    if (A[j] == 1)
			s1 += X[j];
		    else if (A[j] == 2)
			s2 += X[j];
		}
		double p1hat = s1 / (double) n1;
		double p2hat = s2 / (double) n2;

		double rho1hat = 0;
		if (p1hat + p2hat == 0)
		    rho1hat = 0.5;
		else
		    rho1hat = p1hat / (p1hat + p2hat);
		p = 1 - rho1hat;

	    } else
		// NOT SUPPORTED
		p = 0;
	}

	if (p == 0)
	    p = 1/(double)N;
	if (p == 1)
	    p = 1-1/(double)N;

	// DBCD()
	double num1 = (1-p) * pow( (1-p) / (n1/(double)(i+1)) , 2);
	double num2 = p * pow( p / (n2/(double)(i+1)) , 2);
	double phi = num2 / (num1 + num2);

	A[i] = R::rbinom(1, phi) + 1;
	if (A[i] == 1)
	    X[i] = R::rbinom(1, p1);
	else
	    X[i] = R::rbinom(1, p2);
   }
	
   int n0 = 0, n1 = 0, s0 = 0, s1 = 0;
   std::vector<double> x0, x1;
   for (int i = 0; i < A.size(); i++) {
       if (A[i] == 1) {
	   n0++;
	   s0 += X[i];
	   x0.push_back(X[i]);
       } else if (A[i] == 2) { // == 2
	   n1++;
	   s1 += X[i];
	   x1.push_back(X[i]);
       } else {
	   Rcpp::Rcout << "NA BRANCH TAKEN\n";
       }	
    }
    double response = (s0 + s1) / (double) N;
    
    double per_sup;
    if (p1 > p2)
	per_sup = n0 / (double) N;
    else if (p1 < p2)
	per_sup = n1 / (double) N;
    else
	per_sup = NA_REAL;
    
    double Z_P = wald(x0, x1, measure);

    vector<double> out = {response, Z_P, per_sup};
    return out;
	    
}



vector<double> two_arm_ERADE(int N, double p1, double p2, int measure, int burnin, int ar) {
    vector<int> A(N), X(N);

    for (int i = 0; i < burnin; i++) {
	A[i] = 1;
	X[i] = R::rbinom(1, p1);
    }
    for (int i = burnin; i < 2*burnin; i++) {
	A[i] = 2;
	X[i] = R::rbinom(1, p2);
    }
    
    for (int i = 2*burnin; i < N; i++) {
	int n1 = 0, n2 = 0;
	for (int j = 0; j < A.size(); j++) {
	    if (A[j] == 1)
		n1++;
	    if (A[j] == 2)
		n2++;
	}

	double p;
    
	if (ar == NEYMAN) {
	    if (measure == SMD) {
		double sig1_N = 0, sig2_N = 0;
		double mean1 = 0, mean2 = 0;
		int tot1 = 0, tot2 = 0;
		for (int j = 0; j < A.size(); j++) {
		    if (A[j] == 1) {
			mean1 += X[j]; tot1++;
		    } else if (A[j] == 2) {
			mean2 += X[j]; tot2++;
		    }
		}
		mean1 /= tot1; mean2 /= tot2;
		for (int j = 0; j < A.size(); j++) {
		    if (A[j] == 1)
			sig1_N += pow(X[j] - mean1, 2);
		    else if (A[j] == 2)
			sig2_N += pow(X[j] - mean2, 2);
		}
		sig1_N = sqrt(sig1_N / (tot1-1));
		sig2_N = sqrt(sig2_N / (tot2-1));
		
		if (sig1_N + sig2_N == 0)
		    p = 0.5;
		else
		    p = 1 - sig1_N / (sig1_N + sig2_N);
		
	    } else { // LRR
		int s1 = 0, s2 = 0;
		for (int j = 0; j < A.size(); j++) {
		    if (A[j] == 1)
			s1 += X[j];
		    else if (A[j] == 2)
			s2 += X[j];
		}
		double p1hat = s1 / (double) n1;
		double p2hat = s2 / (double) n2;
		if ( (sqrt((1-p1hat)*p2hat) + sqrt((1-p2hat)*p1hat)) == 0)
		    p = 0.5;
		else
		    p = sqrt((1-p1hat)*p2hat) / (sqrt((1-p1hat)*p2hat) + sqrt((1-p2hat)*p1hat));
	    }

	} else if (ar == MINF) {
	    int s1 = 0, s2 = 0;
	    for (int j = 0; j < A.size(); j++) {
		if (A[j] == 1)
		    s1 += X[j];
		else if (A[j] == 2)
		    s2 += X[j];
	    }
	    double p1hat = s1 / (double) n1;
	    double p2hat = s2 / (double) n2;

	    double rho1hat;
	    if (measure == SMD) {
		if (p1hat + p2hat == 0)
		    rho1hat = 0.5;
		else
		    rho1hat = sqrt(p1hat) / (sqrt(p1hat) + sqrt(p2hat));
	    } else { // LRR
		if ((sqrt(p1hat)*(1-p2hat) + sqrt(p2hat)*(1-p1hat)) == 0)
		    rho1hat = 0.5;
		else
		    rho1hat = 1 - sqrt(p2hat)*(1-p1hat) / (sqrt(p1hat)*(1-p2hat) + sqrt(p2hat)*(1-p1hat));
	    }
	    p = 1 - rho1hat;

	} else { // (ar == AD)
	    if (measure == SMD) {
		int s1 = 0, s2 = 0;
		for (int j = 0; j < A.size(); j++) {
		    if (A[j] == 1)
			s1 += X[j];
		    else if (A[j] == 2)
			s2 += X[j];
		}
		double p1hat = s1 / (double) n1;
		double p2hat = s2 / (double) n2;

		double rho1hat = 0;
		if (p1hat + p2hat == 0)
		    rho1hat = 0.5;
		else
		    rho1hat = p1hat / (p1hat + p2hat);
		p = 1 - rho1hat;

	    } else
		// NOT SUPPORTED
		p = 0;
	}

	if (p == 0)
	    p = 1/(double)N;
	if (p == 1)
	    p = 1-1/(double)N;

	double allocprop = n1 / (double)(i+1);

	// ERADE()
	double rhohat = 1-p;
	double phi;
	if (allocprop > rhohat)
	    phi = 0.5 * rhohat;
	else if (allocprop < rhohat)
	    phi = 1 - 0.5 * (1 - rhohat);
	else
	    phi = rhohat;

	A[i] = R::rbinom(1, 1-phi) + 1;
	if (A[i] == 1)
	    X[i] = R::rbinom(1, p1);
	else
	    X[i] = R::rbinom(1, p2);
   }

   int n0 = 0, n1 = 0, s0 = 0, s1 = 0;
   std::vector<double> x0, x1;
   for (int i = 0; i < A.size(); i++) {
       if (A[i] == 1) {
	   n0++;
	   s0 += X[i];
	   x0.push_back(X[i]);
       } else if (A[i] == 2) { // == 2
	   n1++;
	   s1 += X[i];
	   x1.push_back(X[i]);
       } else {
	   Rcpp::Rcout << "NA BRANCH TAKEN\n";
       }	
    }
    double response = (s0 + s1) / (double) N;
    
    double per_sup;
    if (p1 > p2)
	per_sup = n0 / (double) N;
    else if (p1 < p2)
	per_sup = n1 / (double) N;
    else
	per_sup = NA_REAL;
    
    double Z_P = wald(x0, x1, measure);

    vector<double> out = {response, Z_P, per_sup};
    return out;

}

// [[Rcpp::export]]
Rcpp::List simulationCpp(int N, double p1, double p2, int burnin, int nsim, int method, int measure) {
    
    vector<vector<double>> out(nsim);

    switch(method) {
    case ER:
	for (int i = 0; i < nsim; i++)
	    out[i] = two_arm_ER(N, p1, p2, measure, burnin, false);
	break;
    case PW:
	for (int i = 0; i < nsim; i++)
	    out[i] = PWfun(N, p1, p2, burnin, measure, false);
	break;
    case RPW:
	for (int i = 0; i < nsim; i++)
	    out[i] = randomiseRPW(N, p1, p2, burnin, measure, false);
	break;
    case BRAR:
	for (int i = 0; i < nsim; i++)
	    out[i] = BRARfun(N, p1, p2, burnin, measure, false);
	break;
    case SMLE_Neyman:
	for (int i = 0; i < nsim; i++)
	    out[i] = two_arm_SMLE(N, p1, p2, measure, burnin, NEYMAN);
	break;
    case SMLE_minF:
      	for (int i = 0; i < nsim; i++)
	    out[i] = two_arm_SMLE(N, p1, p2, measure, burnin, MINF);
	break;
    case SMLE_AD:
	for (int i = 0; i < nsim; i++)
	    out[i] = two_arm_SMLE(N, p1, p2, measure, burnin, AD);
	break;
    case DBCD_Neyman:
	for (int i = 0; i < nsim; i++)
	    out[i] = two_arm_DBCD(N, p1, p2, measure, burnin, NEYMAN);
	break;
    case DBCD_minF:
	for (int i = 0; i < nsim; i++)
	    out[i] = two_arm_DBCD(N, p1, p2, measure, burnin, MINF);
	break;
    case DBCD_AD:
	for (int i = 0; i < nsim; i++)
	    out[i] = two_arm_DBCD(N, p1, p2, measure, burnin, AD);
	break;
    case ERADE_Neyman:
	for (int i = 0; i < nsim; i++)
	    out[i] = two_arm_ERADE(N, p1, p2, measure, burnin, NEYMAN);
	break;
    case ERADE_minF:
	for (int i = 0; i < nsim; i++)
	    out[i] = two_arm_ERADE(N, p1, p2, measure, burnin, MINF);
	break;
    case ERADE_AD:
	for (int i = 0; i < nsim; i++)
	    out[i] = two_arm_ERADE(N, p1, p2, measure, burnin, AD);
	break;
	
    default:
	Rcpp::Rcout << "Option not supported\n";
    }

    int reject_Z = 0;
    for (int i = 0; i < nsim; i++)
	if (out[i][1] < 0.05)
	    reject_Z++;
    double power = reject_Z / (double) nsim;    
    double power_err = sqrt(power * (1 - power) / nsim);
    power *= 100;
    double per_sup = 0;
    for (int i = 0; i < nsim; i++)
	per_sup += out[i][2];
    per_sup /= nsim;
    double var_per_sup = 0;
    for (int i = 0; i < nsim; i++)
	var_per_sup += pow(out[i][2] - per_sup, 2);
    var_per_sup = var_per_sup / (nsim-1) * 10000;
    per_sup *= 100;
    double emr = 0;
    for (int i = 0; i < nsim; i++)
	emr += out[i][0];
    emr /= nsim;
    double var_emr = 0;
    for (int i = 0; i < nsim; i++)
	var_emr += pow(out[i][0] - emr, 2);
    var_emr /= (nsim-1);
   
    return Rcpp::List::create(power, per_sup, var_per_sup, emr, var_emr, power_err * 100); 

}
    
