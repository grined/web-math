package com.math;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.DecompositionSolver;
import org.apache.commons.math3.linear.LUDecomposition;
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
		for (int i = 0; i <= sol.getMaxIndex(); i++) 
		{
			System.out.print(sol.getEntry(i) + "   ");
		}
	}
	
	public void solve_theory_of_gym()
	{
		double la1 = E1*nu1/(1+nu1)/(1-2*nu1);
		double la2 = E2*nu2/(1+nu2)/(1-2*nu2);
		double la3 = E3*nu3/(1+nu3)/(1-2*nu3);
		double mu1 = E1/2/(1+nu1);
		double mu2 = E2/2/(1+nu2);
		double mu3 = E3/2/(1+nu3);
		double k1 = 2./3*mu1+la1;
		double k2 = 2./3*mu2+la2;
		double k3 = 2./3*mu3+la3;
		double C1 = k*k*E1;
		double C2 = k*k*E2;
		double C3 = k*k*E3;
		
		
		

	}
	
	public void calculate()
	{
		xMax = x2(n);
		xLast = x2(n)+hp;
		solve_system_of_border_cond();
		
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
