import java.io.File;
import java.io.IOException;

// Démonstration du TP2 (parties 1 et 2).
public class Demo {

	static final String DONNEES = "AC_tp2_part2_donnees/";
	static final String RESULTATS = "resultats/";

	public static void main(String[] args) throws IOException {
		new File(RESULTATS).mkdirs();

		partie1();
		partie2_signal1D();
		partie2_images();

		System.out.println("\nFin de la démo. Les images sont dans " + RESULTATS);
	}

	// Partie 1 : FFT, FFT inverse, multiplication de polynômes, temps
	static void partie1() {
		titre("PARTIE 1 : FFT récursive et multiplication de polynômes");

		double[] a = {1, 2, 3, 4};
		CpxTab tfd = FFT_1D.FFT(a);
		System.out.println("TFD(1,2,3,4)                 = " + tfd);
		System.out.println("TFD_inverse(TFD(1,2,3,4))    = " + FFT_1D.FFT_inverse(tfd));

		System.out.println();
		System.out.println("2 x (-3)                     = " + FFT_1D.multiplication_polynome_viaFFT(new double[]{2}, new double[]{-3}));
		System.out.println("(2+X)(-3+2X)                 = " + FFT_1D.multiplication_polynome_viaFFT(new double[]{2, 1}, new double[]{-3, 2}));
		double[] b = {-3, 2, -5, 0};
		System.out.println("(1+2X+3X^2+4X^3)(-3+2X-5X^2) via FFT   = " + FFT_1D.multiplication_polynome_viaFFT(a, b));
		System.out.print(  "(1+2X+3X^2+4X^3)(-3+2X-5X^2) via coeff = ");
		afficher(FFT_1D.multiplication_polynome_viaCoeff(a, b));

		for (int i = 0; i < 20; i++) {
			double[] w1 = FFT_1D.random(1024), w2 = FFT_1D.random(1024);
			FFT_1D.multiplication_polynome_viaCoeff(w1, w2);
			FFT_1D.multiplication_polynome_viaFFT(w1, w2);
		}
		System.out.println("n\tcoeff\tFFT");
		for (int n = 64; n <= 8192; n *= 2) {
			double[] p1 = FFT_1D.random(n), p2 = FFT_1D.random(n);
			long tCoeff = Long.MAX_VALUE, tFFT = Long.MAX_VALUE;
			for (int r = 0; r < 5; r++) {
				long d = System.nanoTime();
				FFT_1D.multiplication_polynome_viaCoeff(p1, p2);
				tCoeff = Math.min(tCoeff, (System.nanoTime() - d) / 1000);

				d = System.nanoTime();
				FFT_1D.multiplication_polynome_viaFFT(p1, p2);
				tFFT = Math.min(tFFT, (System.nanoTime() - d) / 1000);
			}
			System.out.println(n + "\t" + tCoeff + "\t" + tFFT);
		}
	}

	// Partie 2, exercice 1 : TFD de signaux 1D (n = 16)
	static void partie2_signal1D() {
		titre("PARTIE 2, exercice 1 : TFD de signaux 1D (n = 16)");
		int n = 16;
		double[] a = new double[n];

		for (double c : new double[]{3, -2, 0.5}) {
			for (int t = 0; t < n; t++) a[t] = c;
			System.out.println("constante " + c + " : " + FFT_1D.FFT(a));
		}
		for (int k = 0; k < n; k++) {
			for (int t = 0; t < n; t++) a[t] = Math.cos(2 * Math.PI * k * t / n);
			System.out.println("cos k=" + k + " : " + FFT_1D.FFT(a));
		}
		for (int t = 0; t < n; t++) a[t] = Math.cos(2 * Math.PI * t / n) + 0.5 * Math.cos(2 * Math.PI * 3 * t / n);
		System.out.println("cos(1) + 0.5 cos(3) : " + FFT_1D.FFT(a));
		for (int t = 0; t < n; t++) a[t] = 4 + 2 * Math.sin(2 * Math.PI * 2 * t / n) + Math.cos(2 * Math.PI * 7 * t / n);
		System.out.println("4 + 2 sin(2) + cos(7) : " + FFT_1D.FFT(a));
	}

	// Partie 2, exercices 2 à 7 : images
	static void partie2_images() throws IOException {
		CpxImg tigre = charger("tigre_512");
		int n = tigre.taille();

		// Exercice 2 : TFD puis TFD inverse
		titre("PARTIE 2, exercice 2 : TFD puis TFD inverse du tigre");
		CpxImg reconstruit = FFT_2D.FFT_inverse(FFT_2D.FFT(tigre));
		System.out.println("Écart maximal sur les pixels : " + ecartMax(tigre, reconstruit));
		enregistrer(reconstruit, "demo_ex2_tigre_reconstruit");

		// Exercice 3 : spectres
		titre("PARTIE 2, exercice 3 : spectres");
		for (String nom : new String[]{"mire1", "mire2", "mire3", "fingerprint", "barbara_512", "tigre_512"}) {
			CpxImg spectre = FFT_2D.transformeCpxImg(FFT_2D.FFT(charger(nom)));
			enregistrer(spectre, "demo_ex3_spectre_" + nom);
		}

		// Exercice 4 : annulation d'un coefficient
		titre("PARTIE 2, exercice 4 : annulation d'un coefficient du tigre");
		enregistrer(FFT_2D.annuler(tigre, n / 2, n / 2), "demo_ex4_centre");
		enregistrer(FFT_2D.annuler(tigre, n / 2 + 1, n / 2), "demo_ex4_voisin");
		enregistrer(FFT_2D.annuler(tigre, n / 2 + 100, n / 2), "demo_ex4_haute_frequence");

		// Exercice 5 : carré de basses fréquences
		titre("PARTIE 2, exercice 5 : compression par carré de côté 2k (tigre)");
		for (int k = 4; k <= 128; k *= 2) {
			CpxImg F = FFT_2D.FFT(tigre);
			FFT_2D.compression(F, k);
			enregistrer(FFT_2D.FFT_inverse(F), "demo_ex5_k" + k);
			System.out.println("k = " + k + " : " + (4 * k * k) + " coefficients gardés");
		}

		// Exercice 6 : seuillage
		titre("PARTIE 2, exercice 6 : compression par seuillage (tigre, seuil 30)");
		CpxImg F6 = FFT_2D.FFT(tigre);
		int gardes6 = FFT_2D.compression_seuil(F6, 30);
		System.out.println("Coefficients gardés : " + gardes6 + " sur " + (n * n));
		enregistrer(FFT_2D.FFT_inverse(F6), "demo_ex6_seuil30");

		// Exercice 7 : comparaison carré / seuil à 10 %
		titre("PARTIE 2, exercice 7 : carré (k = 81) contre seuil, même nombre de coefficients");
		int cible = 4 * 81 * 81;
		for (String nom : new String[]{"tigre_512", "barbara_512"}) {
			CpxImg I = charger(nom);

			CpxImg Fc = FFT_2D.FFT(I);
			FFT_2D.compression(Fc, 81);
			CpxImg Rc = FFT_2D.FFT_inverse(Fc);

			double seuil = FFT_2D.trouverSeuil(I, cible);
			CpxImg Fs = FFT_2D.FFT(I);
			int gardes = FFT_2D.compression_seuil(Fs, seuil);
			CpxImg Rs = FFT_2D.FFT_inverse(Fs);

			System.out.println(nom + " : carré " + cible + " coefficients, MSE = " + FFT_2D.mse(I, Rc));
			System.out.println(nom + " : seuil " + seuil + ", " + gardes + " coefficients, MSE = " + FFT_2D.mse(I, Rs));
			enregistrer(Rc, "demo_ex7_" + nom + "_carre");
			enregistrer(Rs, "demo_ex7_" + nom + "_seuil");
		}
	}

	static CpxImg charger(String nom) throws IOException {
		return new CpxImg(new BytePixmap(DONNEES + nom + ".pgm"));
	}

	static void enregistrer(CpxImg img, String nom) throws IOException {
		img.convert_to_BytePixmap().write(RESULTATS + nom + ".pgm");
	}

	static double ecartMax(CpxImg a, CpxImg b) {
		double max = 0;
		for (int i = 0; i < a.taille(); i++)
			for (int j = 0; j < a.taille(); j++)
				max = Math.max(max, Math.abs(a.get_p_reel(i, j) - b.get_p_reel(i, j)));
		return max;
	}

	static void afficher(double[] t) {
		StringBuilder s = new StringBuilder("[");
		for (int i = 0; i < t.length; i++) s.append(i > 0 ? ", " : "").append(t[i]);
		System.out.println(s.append("]"));
	}

	static void titre(String t) {
		System.out.println(t);
	}
}
