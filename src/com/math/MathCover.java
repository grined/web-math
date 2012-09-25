package com.math;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.DecompositionSolver;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.QRDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

public class MathCover 
{

/**
 *  Input data
 */
	private double E1,nu1,E2,nu2,E3,nu3,a1,a2,a3,K1,K2,K3,k;
	private double hp,T01,T02;
	private int n;
	private double[] H,h;
//---------------------------------------------------------
	private double xMax, xLast;
	private double[][] A,B,C,F,G;
	
	double la1;
	double la2;
	double la3;
	double mu1;
	double mu2;
	double mu3;
	double k1;
	double k2;
	double k3;
	double C1;
	double C2;
	double C3;

	
	public MathCover()
	{
		H = new double[] {-0, 22,22,22,22,22}; //indexing begin from 1 to n
		h = new double[] {-0, 4,4,4,4,4};
		hp = 2000;
		n = 5;
		E1 = 204;
		nu1 = 0.28;
		E2 = 210;
		nu2 = 0.23;
		E3 = 194;
		nu3 = 0.233;
		a1 = 1.3e-5;
		a2 = 1.03e-5;
		a3 = 1.66e-5;
		K1 = 1./2;
		K2 = 1./90;
		K3 = 1./16;
		k = 1;
		T01 = 2100;
		T02 = 100;
		la1 = E1*nu1/(1+nu1)/(1-2*nu1);
		la2 = E2*nu2/(1+nu2)/(1-2*nu2);
		la3 = E3*nu3/(1+nu3)/(1-2*nu3);
		mu1 = E1/2/(1+nu1);
		mu2 = E2/2/(1+nu2);
		mu3 = E3/2/(1+nu3);
		k1 = 2./3*mu1+la1;
		k2 = 2./3*mu2+la2;
		k3 = 2./3*mu3+la3;
		C1 = k*k*E1;
		C2 = k*k*E2;
		C3 = k*k*E3;
	}
	
	private double x0 (int i)
	{
		double Sum = 0;
		for (int j=1; j<i; j++)
		{
			Sum += (H[j]+h[j]);
		}
		return Sum;
	}

	private double x1 (int i)
	{
		return x0(i)+H[i];
	}
	
	private double x2 (int i)
	{
		return x1(i)+h[i];
	}
	
	
	public void solve_system_of_border_cond() 
	{
		// System of border conditions
		RealMatrix coeff = new Array2DRowRealMatrix(4 * n, 4 * n);
		RealVector constants = new ArrayRealVector(4 * n);
		int constIter = 0;
		// 1st eq
		constants.setEntry(constIter++, T01);
		coeff.addToEntry(0, 0, 1d);
		coeff.addToEntry(0, 1, 0d);
		for (int i = 1; i <= n - 1; i++)
		{
			// 1st
			coeff.addToEntry((i - 1) * 4 + 1, (i - 1) * 2, 1d);
			coeff.addToEntry((i - 1) * 4 + 1, (i - 1) * 2 + 1, x1(i));
			coeff.addToEntry((i - 1) * 4 + 1, 2 * n + (i - 1) * 2, -1d);
			coeff.addToEntry((i - 1) * 4 + 1, 2 * n + (i - 1) * 2 + 1, -x1(i));
			constants.setEntry(constIter++, 0);
			// 2nd
			coeff.addToEntry((i - 1) * 4 + 2, (i - 1) * 2 + 1, K1);
			coeff.addToEntry((i - 1) * 4 + 2, 2 * n + (i - 1) * 2 + 1, -K2);
			constants.setEntry(constIter++, 0);
			// 3rd
			coeff.addToEntry((i - 1) * 4 + 3, (i) * 2, 1d);
			coeff.addToEntry((i - 1) * 4 + 3, (i) * 2 + 1, x2(i));
			coeff.addToEntry((i - 1) * 4 + 3, 2 * n + (i - 1) * 2, -1d);
			coeff.addToEntry((i - 1) * 4 + 3, 2 * n + (i - 1) * 2 + 1, -x2(i));
			constants.setEntry(constIter++, 0);
			// 4th
			coeff.addToEntry((i - 1) * 4 + 4, (i) * 2 + 1, K1);
			coeff.addToEntry((i - 1) * 4 + 4, 2 * n + (i - 1) * 2 + 1, -K2);
			constants.setEntry(constIter++, 0);
		}
		// Last 3 eq
		// 1st
		coeff.addToEntry(4 * n - 3, 2 * n - 2, 1d);
		coeff.addToEntry(4 * n - 3, 2 * n - 1, xMax - h[n]);
		coeff.addToEntry(4 * n - 3, 4 * n - 2, -1d);
		coeff.addToEntry(4 * n - 3, 4 * n - 1, -xMax + h[n]);
		constants.setEntry(constIter++, 0);
		// 2nd
		coeff.addToEntry(4 * n - 2, 2 * n - 1, K1);
		coeff.addToEntry(4 * n - 2, 4 * n - 1, -K2);
		constants.setEntry(constIter++, 0);
		// 3rd
		coeff.addToEntry(4 * n - 1, 4 * n - 2, 1d);
		coeff.addToEntry(4 * n - 1, 4 * n - 1, xMax);
		constants.setEntry(constIter++, 100d);

		DecompositionSolver solver = new LUDecomposition(coeff).getSolver();
		RealVector sol = solver.solve(constants);
		
		C=new double[n+1][];
		G=new double[n+1][];
		int count = sol.getDimension();
		for (int i=1; i<n+1; i++)
		{
			C[i] = new double[3];
			C[i][1] = sol.getEntry(0 + (i-1)*2);
		}
		for (int i=1; i<n+1; i++)
		{
			C[i][2] = sol.getEntry(1 + (i-1)*2);
		}
		for (int i=1; i<n+1; i++)
		{
			G[i] = new double[3];
			G[i][1] = sol.getEntry(2*n + (i-1)*2);
		}
		for (int i=1; i<n+1; i++)
		{
			G[i][2] = sol.getEntry(2*n+1 + (i-1)*2);
		}
	}
	
	private double CE(double C, double E)
	{
		return Math.sqrt(C)/Math.sqrt(E);
	}
	private double E2C(double C, double E)
	{
		return E*E/C;
	}
	
	private double exp(double C, double E, double salt)
	{
		return Math.exp(CE(C,E)*salt);
	}
	
	//Functions
	public double t1(double x, int i)
	{
		return C[i][1] + x*C[i][2];
	}
	
	public double t2(double x, int i)
	{
		return G[i][1] + x*G[i][2];
	}
	
	
	public double r1PART(double x, int i)
	{
		return a1*k1*x*x*C[i][2]/2/E1;
	}
	
	public double r2PART(double x, int i)
	{
		return a2*k2*x*x*G[i][2]/2/E2;
	}
		
	private double A1(double x, int i)
	{
		return exp(C1,E1,x-x1(i))*E1/C1;
	}
	private double A2(double x, int i)
	{
		return exp(C1,E1,x0(i)-x)*E1/C1;
	}
	
	private double A3(double x, int i)
	{
		return 1;
	}
	
	private double A4(double x, int i)
	{
		return x;
	}
	
	private double B1(double x, int i)
	{
		return exp(C2,E2,x-x2(i))*E2/C2;
	}
	
	private double B2(double x, int i)
	{
		return exp(C2,E2,x1(i)-x)*E2/C2;
	}
	
	private double B3(double x, int i)
	{
		return 1;
	}
	
	private double B4(double x, int i)
	{
		return x;
	}
	
	private double F1(double x)
	{
		return exp(C3,E3,x-xLast)*E3/C3;
	}
	private double F2(double x)
	{
		return exp(C3,E3,xMax-x)*E3/C3;
	}
	
	private double F3(double x)
	{
		return 1;
	}
	
	private double F4(double x)
	{
		return x;
	}
	
	// Derivatives 1st level
	public double Dt1(double x, int i)
	{
		return C[i][2];
	}
	
	public double Dt2(double x, int i)
	{
		return G[i][2];
	}
	public double Dr1PART(double x, int i)
	{
		return 2*a1*k1*x*C[i][2]/2/E1;
	}
	
	public double Dr2PART(double x, int i)
	{
		return 2*a2*k2*x*G[i][2]/2/E2;
	}
		
	private double DA1(double x, int i)
	{
		return exp(C1,E1,x-x1(i))*E1/C1*CE(C1,E1);
	}
	
	private double DA2(double x, int i)
	{
		return exp(C1,E1,x0(i)-x)*E1/C1*(-CE(C1,E1));
	}
	
	private double DA3(double x, int i)
	{
		return 0;
	}
	
	private double DA4(double x, int i)
	{
		return 1;
	}
	
	private double DB1(double x, int i)
	{
		return exp(C2,E2,x-x2(i))*E2/C2*CE(C2,E2);
	}
	
	private double DB2(double x, int i)
	{
		return exp(C2,E2,x1(i)-x)*E2/C2*(-CE(C2,E2));
	}
	
	private double DB3(double x, int i)
	{
		return 0;
	}
	
	private double DB4(double x, int i)
	{
		return 1;
	}
	
	private double DF1(double x)
	{
		return exp(C3,E3,x-xLast)*E3/C3*CE(C3,E3);
	}
	private double DF2(double x)
	{
		return exp(C3,E3,xMax-x)*E3/C3*(-CE(C3,E3));
	}
	
	private double DF3(double x)
	{
		return 0;
	}
	
	private double DF4(double x)
	{
		return 1;
	}
	
	//Derivatives 2nd level
	public double DDt1(double x, int i)
	{
		return 0;
	}
	
	public double DDt2(double x, int i)
	{
		return 0;
	}
	public double DDr1PART(double x, int i)
	{
		return 2*a1*k1*C[i][2]/2/E1;
	}
	
	public double DDr2PART(double x, int i)
	{
		return 2*a2*k2*G[i][2]/2/E2;
	}
		
	private double DDA1(double x, int i)
	{
		return exp(C1,E1,x-x1(i))*E1/C1*CE(C1,E1)*CE(C1,E1);
	}
	
	private double DDA2(double x, int i)
	{
		return exp(C1,E1,x0(i)-x)*E1/C1*(-CE(C1,E1))*(-CE(C1,E1));
	}
	
	private double DDA3(double x, int i)
	{
		return 0;
	}
	
	private double DDA4(double x, int i)
	{
		return 0;
	}
	
	private double DDB1(double x, int i)
	{
		return exp(C2,E2,x-x2(i))*E2/C2*CE(C2,E2)*CE(C2,E2);
	}
	
	private double DDB2(double x, int i)
	{
		return exp(C2,E2,x1(i)-x)*E2/C2*(-CE(C2,E2))*(-CE(C2,E2));
	}
	
	private double DDB3(double x, int i)
	{
		return 0;
	}
	
	private double DDB4(double x, int i)
	{
		return 0;
	}
	
	private double DDF1(double x)
	{
		return exp(C3,E3,x-xLast)*E3/C3*CE(C3,E3)*CE(C3,E3);
	}
	private double DDF2(double x)
	{
		return exp(C3,E3,xMax-x)*E3/C3*(-CE(C3,E3))*(-CE(C3,E3));
	}
	
	private double DDF3(double x)
	{
		return 0;
	}
	
	private double DDF4(double x)
	{
		return 0;
	}
	
	//Derivatives 3rd level
	public double DDDt1(double x, int i)
	{
		return 0;
	}
	
	public double DDDt2(double x, int i)
	{
		return 0;
	}
	public double DDDr1PART(double x, int i)
	{
		return 0;
	}
	
	public double DDDr2PART(double x, int i)
	{
		return 0;
	}
		
	private double DDDA1(double x, int i)
	{
		return exp(C1,E1,x-x1(i))*E1/C1*CE(C1,E1)*CE(C1,E1)*CE(C1,E1);
	}
	
	private double DDDA2(double x, int i)
	{
		return exp(C1,E1,x0(i)-x)*E1/C1*(-CE(C1,E1))*(-CE(C1,E1))*(-CE(C1,E1));
	}
	
	private double DDDA3(double x, int i)
	{
		return 0;
	}
	
	private double DDDA4(double x, int i)
	{
		return 0;
	}
	
	private double DDDB1(double x, int i)
	{
		return exp(C2,E2,x-x2(i))*E2/C2*CE(C2,E2)*CE(C2,E2)*CE(C2,E2);
	}
	
	private double DDDB2(double x, int i)
	{
		return exp(C2,E2,x1(i)-x)*E2/C2*(-CE(C2,E2))*(-CE(C2,E2))*(-CE(C2,E2));
	}
	
	private double DDDB3(double x, int i)
	{
		return 0;
	}
	
	private double DDDB4(double x, int i)
	{
		return 0;
	}
	
	private double DDDF1(double x)
	{
		return exp(C3,E3,x-xLast)*E3/C3*CE(C3,E3)*CE(C3,E3)*CE(C3,E3);
	}
	private double DDDF2(double x)
	{
		return exp(C3,E3,xMax-x)*E3/C3*(-CE(C3,E3))*(-CE(C3,E3))*(-CE(C3,E3));
	}
	
	private double DDDF3(double x)
	{
		return 0;
	}
	
	private double DDDF4(double x)
	{
		return 0;
	}
	
	public void solve_theory_of_gym()
	{
		// 2nd system
		RealMatrix coeff = new Array2DRowRealMatrix(8 * n + 5, 8 * n + 4);
		RealVector constants = new ArrayRealVector(8 * n + 5);
		int constIter = 0;
		// 1st eq
		coeff.addToEntry(constIter, 0, DDA1(0, 1));
		coeff.addToEntry(constIter, 1, DDA2(0, 1));
		coeff.addToEntry(constIter, 2, DDA3(0, 1));
		coeff.addToEntry(constIter, 3, DDA4(0, 1));
		constants.setEntry(constIter++, -DDr1PART(0, 1)+a1*Dt1(0, 1));

		// 2nd eq
		coeff.addToEntry(constIter, 0, E1*DA1(0, 1) - E2C(C1,E1)*DDDA1(0,1));
		coeff.addToEntry(constIter, 1, E1*DA2(0, 1) - E2C(C1,E1)*DDDA2(0,1));
		coeff.addToEntry(constIter, 2, E1*DA3(0, 1) - E2C(C1,E1)*DDDA3(0,1));
		coeff.addToEntry(constIter, 3, E1*DA4(0, 1) - E2C(C1,E1)*DDDA4(0,1));
		constants.setEntry(constIter++, -E1*Dr1PART(0, 1)+E2C(C1,E1)*DDDr1PART(0, 1)-a1*(E2C(C1,E1)*DDt1(0, 1)-k1*t1(0, 1)) );
		// 3rd eq
		coeff.addToEntry(constIter, 2, 1d);
		constants.setEntry(constIter++, 0);
		
		double x = 0;
		for (int i = 1; i <= n - 1; i++)
		{
			// 1st
			x = x1(i);
			coeff.addToEntry(constIter, (i - 1) * 4 + 0, A1(x, i));
			coeff.addToEntry(constIter, (i - 1) * 4 + 1, A2(x, i));
			coeff.addToEntry(constIter, (i - 1) * 4 + 2, A3(x, i));
			coeff.addToEntry(constIter, (i - 1) * 4 + 3, A4(x, i));
			
			coeff.addToEntry(constIter, 4 * n + (i - 1) * 4 + 0, -B1(x, i));
			coeff.addToEntry(constIter, 4 * n + (i - 1) * 4 + 1, -B2(x, i));
			coeff.addToEntry(constIter, 4 * n + (i - 1) * 4 + 2, -B3(x, i));
			coeff.addToEntry(constIter, 4 * n + (i - 1) * 4 + 3, -B4(x, i));
			
			constants.setEntry(constIter++, -r1PART(x, i) + r2PART(x, i));			
			// 2nd
			x = x1(i);
			coeff.addToEntry(constIter, (i - 1) * 4 + 0, DA1(x, i));
			coeff.addToEntry(constIter, (i - 1) * 4 + 1, DA2(x, i));
			coeff.addToEntry(constIter, (i - 1) * 4 + 2, DA3(x, i));
			coeff.addToEntry(constIter, (i - 1) * 4 + 3, DA4(x, i));
			
			coeff.addToEntry(constIter, 4 * n + (i - 1) * 4 + 0, -DB1(x, i));
			coeff.addToEntry(constIter, 4 * n + (i - 1) * 4 + 1, -DB2(x, i));
			coeff.addToEntry(constIter, 4 * n + (i - 1) * 4 + 2, -DB3(x, i));
			coeff.addToEntry(constIter, 4 * n + (i - 1) * 4 + 3, -DB4(x, i));
			
			constants.setEntry(constIter++, -Dr1PART(x, i) + Dr2PART(x, i));
			// 3rd
			x = x1(i);
			coeff.addToEntry(constIter, (i - 1) * 4 + 0, DDA1(x, i));
			coeff.addToEntry(constIter, (i - 1) * 4 + 1, DDA2(x, i));
			coeff.addToEntry(constIter, (i - 1) * 4 + 2, DDA3(x, i));
			coeff.addToEntry(constIter, (i - 1) * 4 + 3, DDA4(x, i));
			
			coeff.addToEntry(constIter, 4 * n + (i - 1) * 4 + 0, -DDB1(x, i));
			coeff.addToEntry(constIter, 4 * n + (i - 1) * 4 + 1, -DDB2(x, i));
			coeff.addToEntry(constIter, 4 * n + (i - 1) * 4 + 2, -DDB3(x, i));
			coeff.addToEntry(constIter, 4 * n + (i - 1) * 4 + 3, -DDB4(x, i));
			
			constants.setEntry(constIter++, -DDr1PART(x, i) + DDr2PART(x, i) + a1*Dt1(x, i) - a2*Dt2(x, i));
			//4th
			x = x1(i);
			coeff.addToEntry(constIter, (i - 1) * 4 + 0, E1*DA1(x, i)-E2C(E1,C1)*DDDA1(x, i));
			coeff.addToEntry(constIter, (i - 1) * 4 + 1, E1*DA2(x, i)-E2C(E1,C1)*DDDA2(x, i));
			coeff.addToEntry(constIter, (i - 1) * 4 + 2, E1*DA3(x, i)-E2C(E1,C1)*DDDA3(x, i));
			coeff.addToEntry(constIter, (i - 1) * 4 + 3, E1*DA4(x, i)-E2C(E1,C1)*DDDA4(x, i));
			
			coeff.addToEntry(constIter, 4 * n + (i - 1) * 4 + 0, -E2*DB1(x, i)+E2C(E2,C2)*DDDB1(x, i));
			coeff.addToEntry(constIter, 4 * n + (i - 1) * 4 + 1, -E2*DB2(x, i)+E2C(E2,C2)*DDDB2(x, i));
			coeff.addToEntry(constIter, 4 * n + (i - 1) * 4 + 2, -E2*DB3(x, i)+E2C(E2,C2)*DDDB3(x, i));
			coeff.addToEntry(constIter, 4 * n + (i - 1) * 4 + 3, -E2*DB4(x, i)+E2C(E2,C2)*DDDB4(x, i));
			
			constants.setEntry(constIter++, -E1*Dr1PART(x, i)+E2C(E1,C1)*DDDr1PART(x, i)+E2*Dr2PART(x, i)-E2C(E2,C2)*DDDr2PART(x, i)
											-a1*(E2C(E1,C1)*DDt1(x, i)-k1*t1(x, i))
											+a2*(E2C(E2,C2)*DDt2(x, i)-k2*t2(x, i)));
			// 5th
			x = x2(i);
			coeff.addToEntry(constIter, (i) * 4 + 0, A1(x, i+1));
			coeff.addToEntry(constIter, (i) * 4 + 1, A2(x, i+1));
			coeff.addToEntry(constIter, (i) * 4 + 2, A3(x, i+1));
			coeff.addToEntry(constIter, (i) * 4 + 3, A4(x, i+1));
			
			coeff.addToEntry(constIter, 4 * n + (i - 1) * 4 + 0, -B1(x, i));
			coeff.addToEntry(constIter, 4 * n + (i - 1) * 4 + 1, -B2(x, i));
			coeff.addToEntry(constIter, 4 * n + (i - 1) * 4 + 2, -B3(x, i));
			coeff.addToEntry(constIter, 4 * n + (i - 1) * 4 + 3, -B4(x, i));
			
			constants.setEntry(constIter++, -r1PART(x, i+1) + r2PART(x, i));			
			// 6th
			x = x2(i);
			coeff.addToEntry(constIter, (i) * 4 + 0, DA1(x, i+1));
			coeff.addToEntry(constIter, (i) * 4 + 1, DA2(x, i+1));
			coeff.addToEntry(constIter, (i) * 4 + 2, DA3(x, i+1));
			coeff.addToEntry(constIter, (i) * 4 + 3, DA4(x, i+1));
			
			coeff.addToEntry(constIter, 4 * n + (i - 1) * 4 + 0, -DB1(x, i));
			coeff.addToEntry(constIter, 4 * n + (i - 1) * 4 + 1, -DB2(x, i));
			coeff.addToEntry(constIter, 4 * n + (i - 1) * 4 + 2, -DB3(x, i));
			coeff.addToEntry(constIter, 4 * n + (i - 1) * 4 + 3, -DB4(x, i));
			
			constants.setEntry(constIter++, -Dr1PART(x, i+1) + Dr2PART(x, i));
			// 7th
			x = x2(i);
			coeff.addToEntry(constIter, (i) * 4 + 0, DDA1(x, i+1));
			coeff.addToEntry(constIter, (i) * 4 + 1, DDA2(x, i+1));
			coeff.addToEntry(constIter, (i) * 4 + 2, DDA3(x, i+1));
			coeff.addToEntry(constIter, (i) * 4 + 3, DDA4(x, i+1));
			
			coeff.addToEntry(constIter, 4 * n + (i - 1) * 4 + 0, -DDB1(x, i));
			coeff.addToEntry(constIter, 4 * n + (i - 1) * 4 + 1, -DDB2(x, i));
			coeff.addToEntry(constIter, 4 * n + (i - 1) * 4 + 2, -DDB3(x, i));
			coeff.addToEntry(constIter, 4 * n + (i - 1) * 4 + 3, -DDB4(x, i));
			
			constants.setEntry(constIter++, -DDr1PART(x, i+1) + DDr2PART(x, i) + a1*Dt1(x, i) - a2*Dt2(x, i));
			//8th
			x = x2(i);
			coeff.addToEntry(constIter, (i) * 4 + 0, E1*DA1(x, i+1)-E2C(E1,C1)*DDDA1(x, i+1));
			coeff.addToEntry(constIter, (i) * 4 + 1, E1*DA2(x, i+1)-E2C(E1,C1)*DDDA2(x, i+1));
			coeff.addToEntry(constIter, (i) * 4 + 2, E1*DA3(x, i+1)-E2C(E1,C1)*DDDA3(x, i+1));
			coeff.addToEntry(constIter, (i) * 4 + 3, E1*DA4(x, i+1)-E2C(E1,C1)*DDDA4(x, i+1));
			
			coeff.addToEntry(constIter, 4 * n + (i - 1) * 4 + 0, -E2*DB1(x, i)+E2C(E2,C2)*DDDB1(x, i));
			coeff.addToEntry(constIter, 4 * n + (i - 1) * 4 + 1, -E2*DB2(x, i)+E2C(E2,C2)*DDDB2(x, i));
			coeff.addToEntry(constIter, 4 * n + (i - 1) * 4 + 2, -E2*DB3(x, i)+E2C(E2,C2)*DDDB3(x, i));
			coeff.addToEntry(constIter, 4 * n + (i - 1) * 4 + 3, -E2*DB4(x, i)+E2C(E2,C2)*DDDB4(x, i));
			
			constants.setEntry(constIter++, -E1*Dr1PART(x, i+1)+E2C(E1,C1)*DDDr1PART(x, i+1)+E2*Dr2PART(x, i)-E2C(E2,C2)*DDDr2PART(x, i)
											-a1*(E2C(E1,C1)*DDt1(x, i+1)-k1*t1(x, i+1))
											+a2*(E2C(E2,C2)*DDt2(x, i)-k2*t2(x, i)));
		}
		
		// Last 10 equestions
		int i = n;
		// 1st
		x = xMax-h[n];
		coeff.addToEntry(constIter, (i - 1) * 4 + 0, A1(x, i));
		coeff.addToEntry(constIter, (i - 1) * 4 + 1, A2(x, i));
		coeff.addToEntry(constIter, (i - 1) * 4 + 2, A3(x, i));
		coeff.addToEntry(constIter, (i - 1) * 4 + 3, A4(x, i));
		
		coeff.addToEntry(constIter, 4 * n + (i - 1) * 4 + 0, -B1(x, i));
		coeff.addToEntry(constIter, 4 * n + (i - 1) * 4 + 1, -B2(x, i));
		coeff.addToEntry(constIter, 4 * n + (i - 1) * 4 + 2, -B3(x, i));
		coeff.addToEntry(constIter, 4 * n + (i - 1) * 4 + 3, -B4(x, i));
		
		constants.setEntry(constIter++, -r1PART(x, i) + r2PART(x, i));		
		// 2nd
		x = xMax-h[n];
		coeff.addToEntry(constIter, (i - 1) * 4 + 0, DA1(x, i));
		coeff.addToEntry(constIter, (i - 1) * 4 + 1, DA2(x, i));
		coeff.addToEntry(constIter, (i - 1) * 4 + 2, DA3(x, i));
		coeff.addToEntry(constIter, (i - 1) * 4 + 3, DA4(x, i));
		
		coeff.addToEntry(constIter, 4 * n + (i - 1) * 4 + 0, -DB1(x, i));
		coeff.addToEntry(constIter, 4 * n + (i - 1) * 4 + 1, -DB2(x, i));
		coeff.addToEntry(constIter, 4 * n + (i - 1) * 4 + 2, -DB3(x, i));
		coeff.addToEntry(constIter, 4 * n + (i - 1) * 4 + 3, -DB4(x, i));
		
		constants.setEntry(constIter++, -Dr1PART(x, i) + Dr2PART(x, i));
		// 3rd
		x = xMax-h[n];
		coeff.addToEntry(constIter, (i - 1) * 4 + 0, DDA1(x, i));
		coeff.addToEntry(constIter, (i - 1) * 4 + 1, DDA2(x, i));
		coeff.addToEntry(constIter, (i - 1) * 4 + 2, DDA3(x, i));
		coeff.addToEntry(constIter, (i - 1) * 4 + 3, DDA4(x, i));
		
		coeff.addToEntry(constIter, 4 * n + (i - 1) * 4 + 0, -DDB1(x, i));
		coeff.addToEntry(constIter, 4 * n + (i - 1) * 4 + 1, -DDB2(x, i));
		coeff.addToEntry(constIter, 4 * n + (i - 1) * 4 + 2, -DDB3(x, i));
		coeff.addToEntry(constIter, 4 * n + (i - 1) * 4 + 3, -DDB4(x, i));
		
		constants.setEntry(constIter++, -DDr1PART(x, i) + DDr2PART(x, i) + a1*Dt1(x, i) - a2*Dt2(x, i));
		//4th
		x = xMax-h[n];
		coeff.addToEntry(constIter, (i - 1) * 4 + 0, E1*DA1(x, i)-E2C(E1,C1)*DDDA1(x, i));
		coeff.addToEntry(constIter, (i - 1) * 4 + 1, E1*DA2(x, i)-E2C(E1,C1)*DDDA2(x, i));
		coeff.addToEntry(constIter, (i - 1) * 4 + 2, E1*DA3(x, i)-E2C(E1,C1)*DDDA3(x, i));
		coeff.addToEntry(constIter, (i - 1) * 4 + 3, E1*DA4(x, i)-E2C(E1,C1)*DDDA4(x, i));
		
		coeff.addToEntry(constIter, 4 * n + (i - 1) * 4 + 0, -E2*DB1(x, i)+E2C(E2,C2)*DDDB1(x, i));
		coeff.addToEntry(constIter, 4 * n + (i - 1) * 4 + 1, -E2*DB2(x, i)+E2C(E2,C2)*DDDB2(x, i));
		coeff.addToEntry(constIter, 4 * n + (i - 1) * 4 + 2, -E2*DB3(x, i)+E2C(E2,C2)*DDDB3(x, i));
		coeff.addToEntry(constIter, 4 * n + (i - 1) * 4 + 3, -E2*DB4(x, i)+E2C(E2,C2)*DDDB4(x, i));
		
		constants.setEntry(constIter++, -E1*Dr1PART(x, i)+E2C(E1,C1)*DDDr1PART(x, i)+E2*Dr2PART(x, i)-E2C(E2,C2)*DDDr2PART(x, i)
										-a1*(E2C(E1,C1)*DDt1(x, i)-k1*t1(x, i))
										+a2*(E2C(E2,C2)*DDt2(x, i)-k2*t2(x, i)));
		// 5th
		x = xMax;		
		coeff.addToEntry(constIter, 4 * n + (i - 1) * 4 + 0, B1(x, i));
		coeff.addToEntry(constIter, 4 * n + (i - 1) * 4 + 1, B2(x, i));
		coeff.addToEntry(constIter, 4 * n + (i - 1) * 4 + 2, B3(x, i));
		coeff.addToEntry(constIter, 4 * n + (i - 1) * 4 + 3, B4(x, i));
		
		coeff.addToEntry(constIter, 8 * n + 0, -F1(x));
		coeff.addToEntry(constIter, 8 * n + 1, -F2(x));
		coeff.addToEntry(constIter, 8 * n + 2, -F3(x));
		coeff.addToEntry(constIter, 8 * n + 3, -F4(x));

		constants.setEntry(constIter++, -r2PART(x, i));		
		// 6th
		x = xMax;
		coeff.addToEntry(constIter, 4 * n + (i - 1) * 4 + 0, DB1(x, i));
		coeff.addToEntry(constIter, 4 * n + (i - 1) * 4 + 1, DB2(x, i));
		coeff.addToEntry(constIter, 4 * n + (i - 1) * 4 + 2, DB3(x, i));
		coeff.addToEntry(constIter, 4 * n + (i - 1) * 4 + 3, DB4(x, i));

		coeff.addToEntry(constIter, 8 * n + 0, -DF1(x));
		coeff.addToEntry(constIter, 8 * n + 1, -DF2(x));
		coeff.addToEntry(constIter, 8 * n + 2, -DF3(x));
		coeff.addToEntry(constIter, 8 * n + 3, -DF4(x));
		
		
		constants.setEntry(constIter++, -Dr2PART(x, i));
		// 7th
		x = xMax;
		coeff.addToEntry(constIter, 4 * n + (i - 1) * 4 + 0, DDB1(x, i));
		coeff.addToEntry(constIter, 4 * n + (i - 1) * 4 + 1, DDB2(x, i));
		coeff.addToEntry(constIter, 4 * n + (i - 1) * 4 + 2, DDB3(x, i));
		coeff.addToEntry(constIter, 4 * n + (i - 1) * 4 + 3, DDB4(x, i));

		coeff.addToEntry(constIter, 8 * n + 0, -DDF1(x));
		coeff.addToEntry(constIter, 8 * n + 1, -DDF2(x));
		coeff.addToEntry(constIter, 8 * n + 2, -DDF3(x));
		coeff.addToEntry(constIter, 8 * n + 3, -DDF4(x));
		
		
		constants.setEntry(constIter++, -DDr2PART(x, i) + a2*Dt2(x, i));
		//8th
		x = xMax;
		coeff.addToEntry(constIter, 4 * n + (i - 1) * 4 + 0, E2*DB1(x, i)-E2C(E2,C2)*DDDB1(x, i));
		coeff.addToEntry(constIter, 4 * n + (i - 1) * 4 + 1, E2*DB2(x, i)-E2C(E2,C2)*DDDB2(x, i));
		coeff.addToEntry(constIter, 4 * n + (i - 1) * 4 + 2, E2*DB3(x, i)-E2C(E2,C2)*DDDB3(x, i));
		coeff.addToEntry(constIter, 4 * n + (i - 1) * 4 + 3, E2*DB4(x, i)-E2C(E2,C2)*DDDB4(x, i));

		coeff.addToEntry(constIter, 8 * n + 0, -E3*DF1(x) + E2C(E3,C3)*DDDF1(x));
		coeff.addToEntry(constIter, 8 * n + 1, -E3*DF2(x) + E2C(E3,C3)*DDDF2(x));
		coeff.addToEntry(constIter, 8 * n + 2, -E3*DF3(x) + E2C(E3,C3)*DDDF3(x));
		coeff.addToEntry(constIter, 8 * n + 3, -E3*DF4(x) + E2C(E3,C3)*DDDF4(x));
	
		constants.setEntry(constIter++, -E2*Dr2PART(x, i)+E2C(E2,C2)*DDDr2PART(x, i)
										-a2*(E2C(E2,C2)*DDt2(x, i)-k2*t2(x, i))
										-a3*k3*100);
		//9th
		x = hp + xMax;

		coeff.addToEntry(constIter, 8 * n + 0, DDF1(x));
		coeff.addToEntry(constIter, 8 * n + 1, DDF2(x));
		coeff.addToEntry(constIter, 8 * n + 2, DDF3(x));
		coeff.addToEntry(constIter, 8 * n + 3, DDF4(x));
	
		constants.setEntry(constIter++, 0);

		//10th
		x = hp + xMax;

		coeff.addToEntry(constIter, 8 * n + 0, E3*DF1(x) - E2C(E3,C3)*DDDF1(x));
		coeff.addToEntry(constIter, 8 * n + 1, E3*DF2(x) - E2C(E3,C3)*DDDF2(x));
		coeff.addToEntry(constIter, 8 * n + 2, E3*DF3(x) - E2C(E3,C3)*DDDF3(x));
		coeff.addToEntry(constIter, 8 * n + 3, E3*DF4(x) - E2C(E3,C3)*DDDF4(x));
	
		constants.setEntry(constIter++, a3*k3*100);

		DecompositionSolver solver = new QRDecomposition(coeff).getSolver();
		RealVector sol = solver.solve(constants);
//		for (int j=0; j<sol.getDimension(); j++)	
//			System.out.print(sol.getEntry(j)+ "  ");
		for (int j=0; j<8*n+5; j++)
		{
			for (int l=0; l<8*n+4; l++)
			{
				if ((Math.abs(coeff.getEntry(j,l)-MathCoverTest.test[j][l]))<0.0001)
				{
					System.out.print(j+"  "+l);
					return;
				}
			}
		}
	}
	
	public void calculate()
	{
		xMax = x2(n);
		xLast = x2(n)+hp;
		solve_system_of_border_cond();
		solve_theory_of_gym();
		
	}
	
	public double getE1() {
		return E1;
	}

	public void setE1(double e1) {
		E1 = e1;
	}

	public double getNu1() {
		return nu1;
	}

	public void setNu1(double nu1) {
		this.nu1 = nu1;
	}

	public double getE2() {
		return E2;
	}

	public void setE2(double e2) {
		E2 = e2;
	}

	public double getNu2() {
		return nu2;
	}

	public void setNu2(double nu2) {
		this.nu2 = nu2;
	}

	public double getE3() {
		return E3;
	}

	public void setE3(double e3) {
		E3 = e3;
	}

	public double getNu3() {
		return nu3;
	}

	public void setNu3(double nu3) {
		this.nu3 = nu3;
	}

	public double getA1() {
		return a1;
	}

	public void setA1(double a1) {
		this.a1 = a1;
	}

	public double getA2() {
		return a2;
	}

	public void setA2(double a2) {
		this.a2 = a2;
	}

	public double getA3() {
		return a3;
	}

	public void setA3(double a3) {
		this.a3 = a3;
	}

	public double getK1() {
		return K1;
	}

	public void setK1(double k1) {
		K1 = k1;
	}

	public double getK2() {
		return K2;
	}

	public void setK2(double k2) {
		K2 = k2;
	}

	public double getK() {
		return k;
	}

	public void setK(double k) {
		this.k = k;
	}

	public double getHp() {
		return hp;
	}

	public void setHp(double hp) {
		this.hp = hp;
	}

	public int getN() {
		return n;
	}

	public void setN(int n) {
		this.n = n;
	}

	public double getT01() {
		return T01;
	}

	public void setT01(double t01) {
		T01 = t01;
	}

	public double getT02() {
		return T02;
	}

	public void setT02(double t02) {
		T02 = t02;
	}
	
	public static void main(String[] argv)
	{
		System.out.println();
		new MathCover().calculate();
	}
}
