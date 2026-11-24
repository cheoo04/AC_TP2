import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Locale;
import java.util.Random;
import java.util.Scanner;
/**** on va ici implémenter la transformée de Fourier rapide 1D ****/

public class FFT_1D {
	
	//"combine" c1 et c2 selon la formule vue en TD
	// c1 et c2 sont de même taille
	// la taille du résultat est le double de la taille de c1
	public static CpxTab combine(CpxTab c1, CpxTab c2) {
		assert (c1.taille()==c2.taille()) : "combine: c1 et c2 ne sont pas de même taille, taille c1="+c1.taille()+" taille c2="+c2.taille();

		int m = c1.taille();
		
		int n = 2*m;

		CpxTab res = new CpxTab(n);

		for (int k = 0; k<m ; k++) {
			
			double theta = 2.0 * Math.PI * k / n;
			double wr = Math.cos(theta);
			double wi = Math.sin(theta);
			
			double r_c1 = c1.get_p_reel(k);
			double i_c1 = c1.get_p_imag(k);

			double r_c2 = c2.get_p_reel(k);
			double i_c2 = c2.get_p_imag(k);

			double r_prod = wr*r_c2 - wi*i_c2;
			double i_prod = wr*i_c2 + wi*r_c2;

			res.set_p_reel(k, r_c1 + r_prod);
			res.set_p_imag(k, i_c1 + i_prod);
			res.set_p_reel(k + m, r_c1 - r_prod);
			res.set_p_imag(k + m, i_c1 - i_prod);
		}


		return res;
	}

	//renvoie la TFD d'un tableau de complexes
	//la taille de x doit être une puissance de 2
	public static CpxTab FFT(CpxTab x) {
		
		if ( x.taille()==1 ) { return new CpxTab(x); }

		assert (x.taille()%2==0) : "FFT: la taille de x doit être une puissance de 2";
		
		int m = x.taille() / 2;

		CpxTab pair = new CpxTab(m);
		CpxTab impair = new CpxTab(m);

		for (int k=0; k<m; k++){
			pair.set_p_reel(k, x.get_p_reel(2*k));
			pair.set_p_imag(k, x.get_p_imag(2*k));

			impair.set_p_reel(k, x.get_p_reel(2*k + 1));
			impair.set_p_imag(k, x.get_p_imag(2*k + 1));
		}

		CpxTab c1 = FFT(pair);
		CpxTab c2 = FFT(impair);

		return combine(c1, c2);
	}

	//renvoie la TFD d'un tableau de réels
	//la taille de x doit être une puissance de 2
	public static CpxTab FFT(double[] x) {
		return FFT(new CpxTab(x));
	}
			
	//renvoie la transformée de Fourier inverse de y
	public static CpxTab FFT_inverse(CpxTab y) {
		
		int n = y.taille();

		CpxTab conjug = y.conjugue();
		CpxTab fft_y = FFT(conjug);
		CpxTab res = fft_y.conjugue();

		for (int k = 0; k<n; k++){
			res.set_p_reel(k, res.get_p_reel(k) / n);
			res.set_p_imag(k, res.get_p_imag(k) / n);
		}

		return res;
	}
	
	//calcule le produit de deux polynômes en utilisant la FFT
	//tab1 et tab2, sont les coefficients de ces polynômes
	// CpxTab sera le tableau des coefficients du polynôme produit (purement réel)
	public static CpxTab multiplication_polynome_viaFFT(double[] tab1, double[] tab2) {
		
		//on commence par doubler la taille des vecteurs en rajoutant des zéros à la fin (cf TD)
		double[] t1 = new double[2*tab1.length], t2 = new double[2*tab1.length];
		for (int i = 0; i < tab1.length; i++) {
			t1[i] = tab1[i];
			t2[i] = tab2[i];
		}

		CpxTab c1 = FFT(t1);
		CpxTab c2 = FFT(t2);

		CpxTab prod = CpxTab.multiplie(c1, c2);

		return FFT_inverse(prod);
	}

	
	//renvoie un tableau de réels aléatoires
	//utile pour tester la multiplication de polynômes
	public static double[] random(int n) {
		double[] t = new double[n];

		for (int i = 0; i < n; i++)
			t[i] = Math.random();
		return t;
	}

	//effectue la multiplication de polynômes représentés par coefficients
	// p1, p2 les coefficients des deux polynômes P1 et P2
	// renvoie les coefficients du polynôme P1*P2
	public static double [] multiplication_polynome_viaCoeff(double[] p1, double[] p2){
		
		int n = p1.length + p2.length - 1;
		double a,b;
		double [] out = new double[n];
		for (int k = 0; k < n; k++) {
			for (int i = 0; i <= k; i++) {
				a = (i<p1.length) ? p1[i]:0;
				b = (k-i<p2.length) ? p2[k-i] : 0;
				out[k] += a*b;
			}
		}
		return out;
	}
	

	//affiche un tableau de réels
	private static void afficher(double [] t){
		System.out.print("[");
		for(int k=0;k<t.length;k++){
			System.out.print(t[k]);
			if (k<(t.length-1))
				System.out.print(" ");
		}
		System.out.println("]");
	}

	private static void afficherTFD(String nom, double[] a) {
    System.out.println(nom + " : " + FFT(a));
	}

	
	public static void main(String[] args) {
		// double[] t5 = {1,2,3,4};

		// /* Exo 2: calculez et affichez TFD(1,2,3,4) */
		// 	CpxTab TFD_t5 = FFT(t5);
		// 	System.out.println(" TFD(1,2,3,4) = " + TFD_t5);
		
		// /* Exo 3: calculez et affichez TFD_inverse(TFD(1,2,3,4)) */
		// 	CpxTab TFD_inverse_t5 = FFT_inverse( TFD_t5);
		// 	System.out.println(" TFD_inverse(TFD(1,2,3,4)) = " + TFD_inverse_t5);		

		// /* Exo 4: multiplication polynomiale, vérification*/
		// 	/* A(X) = 2 et B(X)=-3 */
		// 	double[] t1 = {2};
		// 	double[] t2 = {-3};
		// 	CpxTab mult1 = multiplication_polynome_viaFFT(t1, t2); 
		// 	System.out.println("La multiplication polynomiale de A(X) = 2 et B(X)=-3 est : C(X)= " + mult1);


		// 	/* A(X) = 2+X et B(X)= -3+2X */
		// 	double[] t3 = {2, 1};
		// 	double[] t4 = {-3, 2};
		// 	CpxTab mult2 = multiplication_polynome_viaFFT(t3, t4); 
		// 	System.out.println("La multiplication polynomiale de A(X) = 2+X et B(X)= -3+2X est : C(X)= " + mult2);
					

		// 	/* A(X) = 1 + 2X + 3X^2 + 4X^3 et B(X) = -3 + 2X - 5 X^2*/
		// 	System.out.println("-----------------------------------------------------");
		// 	System.out.println("   Comparaison des 2 méthodes de multiplications polynomiales");
		// 	double[] t6 = {-3,2,-5,0};
		// 	System.out.println("mult via FFT  --> " + multiplication_polynome_viaFFT(t5, t6));
		// 	System.out.print(  "mult via coeff -> ");
		// 	afficher(multiplication_polynome_viaCoeff(t5, t6));

		// /* Exo 5: comparaison des temps de calculs */
	
		// // Pour étude du temps de calcul
		// // Echauffement 
		// 	for (int i = 0; i < 20; i++) {
		// 		double[] w1 = random(1024), w2 = random(1024);
		// 		multiplication_polynome_viaCoeff(w1, w2);
		// 		multiplication_polynome_viaFFT(w1, w2);
		// 	}

		// // Coeff
		// 	double a=0; 
		// 	List<Integer> tailles_Coeff = new ArrayList<Integer>();
		// 	List<Integer> temps_Coeff = new ArrayList<Integer>();
		// 	for (int n = 64; n<= 32768; n = n*2){  // taille des polynômes à multiplier (testez différentes valeurs en gardant des puissances de 2)
			
		// 		tailles_Coeff.add(n);
		// 		double[] tab1 =random(n),tab2 = random(n);

		// 		long meilleur = Long.MAX_VALUE;
		// 		for (int r = 0; r < 5; r++) {
		// 			long date1 = System.nanoTime();
		// 			multiplication_polynome_viaCoeff(tab1, tab2);
		// 			long date2 = System.nanoTime();
		// 			long date = (date2 - date1) / 1000;
		// 			meilleur = Math.min(meilleur, date);
		// 		}
		// 		temps_Coeff.add((int) meilleur);

		// 	}
		// 	System.out.println("   via Coeff: " + temps_Coeff);
			
		// 	System.out.println("n \t T(n) \t T/n^2 \t T(n)/T(n/2)");
		// 	for (int i = 0; i < tailles_Coeff.size(); i++) {
		// 		double n = tailles_Coeff.get(i);
		// 		double t = temps_Coeff.get(i);

		// 		double rapportModele = t / (n * n);

		// 		String rapportDouble = "-";
		// 		if (i > 0) {
		// 			double tPrecedent = temps_Coeff.get(i - 1);
		// 			rapportDouble = String.valueOf(t / tPrecedent);
		// 		}

		// 		System.out.println(n + "\t" + t + "\t" + rapportModele + "\t" + rapportDouble);
		// 	}


		// 	Graph g1 = new Graph(tailles_Coeff, temps_Coeff);
		// 	g1.display();
		// 	Scanner in = new Scanner(System.in);
		// 	in.useLocale(Locale.US);
		// 	System.out.print("a = ");
		// 	a = in.nextDouble();
		// 	g1.addQuadratic(a);

		// // FFT
		// 	List<Integer> tailles_FFT = new ArrayList<Integer>();
		// 	List<Integer> temps_FFT = new ArrayList<Integer>();
		// 	for (int n = 64; n<= 32768; n = n*2){  // taille des polynômes à multiplier (testez différentes valeurs en gardant des puissances de 2)
			
		// 		tailles_FFT.add(n);
		// 		double[] tab1 =random(n),tab2 = random(n);
		// 		long meilleur = Long.MAX_VALUE;
		// 		for (int r = 0; r < 5; r++) {
		// 			long date1 = System.nanoTime();
		// 			multiplication_polynome_viaFFT(tab1, tab2);
		// 			long date2 = System.nanoTime();
		// 			long date = (date2 - date1) / 1000;
		// 			meilleur = Math.min(meilleur, date);
		// 		}
		// 		temps_FFT.add((int) meilleur);

		// 	}
		// 	System.out.println("   via FFT  : " + temps_FFT);

		// 	System.out.println("n \t T(n) \t T/nln(n) \t T(n)/T(n/2)");
		// 	for (int i = 0; i < tailles_FFT.size(); i++) {
		// 		double n = tailles_FFT.get(i);
		// 		double t = temps_FFT.get(i);

		// 		double rapportModele = t / (n * Math.log(n));

		// 		String rapportDouble = "-";
		// 		if (i > 0) {
		// 			double tPrecedent = temps_FFT.get(i - 1);
		// 			rapportDouble = String.valueOf(t / tPrecedent);
		// 		}

		// 		System.out.println(n + "\t" + t + "\t" + rapportModele + "\t" + rapportDouble);
		// 	}

		// 	Graph g2 = new Graph(tailles_FFT, temps_FFT);
		// 	g2.display();
		// 	System.out.print("a = ");
		// 	a = in.nextDouble();
		// 	g2.addnLogn(a);
		// 	in.close();

		int n = 16;
		double[] a = new double[n];

		// 1) vecteur constant, pour plusieurs valeurs de la constante
		double[] constantes = {3, -2, 0.5};
		for (double c : constantes) {
			for (int t = 0; t < n; t++) {
				a[t] = c;
			}
			afficherTFD("constante a=" + c, a);
		}

		// 2) sinusoïde pure, pour toutes les fréquences k
		for (int k = 0; k < n; k++) {
			for (int t = 0; t < n; t++) {
				a[t] = Math.cos(2 * Math.PI * k * t / n);
			}
			afficherTFD("cos k=" + k, a);
		}

		// 3) somme de 2 sinusoïdes
		for (int t = 0; t < n; t++) {
			a[t] = Math.cos(2 * Math.PI * t / n) + 0.5 * Math.cos(2 * Math.PI * 3 * t / n);
		}
		afficherTFD("cos(1) + 0.5 cos(3)", a);

		// 4) la fonction 4 + 2 sin + cos
		for (int t = 0; t < n; t++){
			a[t] = 4 + 2 * Math.sin(2 * Math.PI * 2 * t / n) + Math.cos(2 * Math.PI * 7 * t / n);
		}
		afficherTFD("4 + 2 sin(2) + cos(7)", a);


	}

}