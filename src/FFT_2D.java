import java.io.IOException;

public class FFT_2D {

	//renvoie la TFD d'une image de complexes
	public static CpxImg FFT(CpxImg I) {

		CpxImg out = new CpxImg(I.taille());

		// FFT 1D sur les lignes
		for (int k = 0; k < I.taille(); k++)
			out.set_line(k,FFT_1D.FFT(I.get_line(k)));
		  
		// transposition
		out.transpose();

		// FFT 1D sur les "nouvelles" lignes de out (anciennes colonnes)
		for (int k = 0; k < I.taille(); k++)
			out.set_line(k,FFT_1D.FFT(out.get_line(k)));

		//on re transpose pour revenir dans le sens de départ
		out.transpose();
		
		//on divise par la taille de I
		out.multiply(1./I.taille());
		return out.recentrage();
	}
	
	//renvoie la TFD inverse d'une images de complexes
	public static CpxImg FFT_inverse(CpxImg I) {
		I = I.recentrage();
		CpxImg out = new CpxImg(I.taille());
		for (int k = 0; k < I.taille(); k++)
			out.set_line(k, I.get_line(k).conjugue());

		out = FFT(out).recentrage();
		for (int k = 0; k < I.taille(); k++)
			out.set_line(k, out.get_line(k).conjugue());
		return out;
	}

	// compression par mise à zéro des coefficients de fréquence 
	// FI contient la TDF de I 
	// Dans FI on met à zéros tous les coefficients correspondant à des fréquences inférieures à k
	public static void compression(CpxImg FI, int k) {
		int n = FI.taille();
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				if (i <= n/2 - k || i > n/2 + k || j <= n/2 - k || j > n/2 + k  ){
					FI.set_p_reel(i, j, 0);
					FI.set_p_imag(i, j, 0);
				}
			}
		}
	}

	// compression par seuillage des coefficients faibles
	// FI contient la TDF de I 
	// Dans FI on met à zéros tous les coefficients dont le module est inférieur à seuil 
	// on renvoie le nombre de coefficients conservés 
	public static int compression_seuil(CpxImg FI, double seuil){
		
		int n = FI.taille();
		int gardes = 0;
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				double re = FI.get_p_reel(i, j);
				double im = FI.get_p_imag(i, j);
				double mod = Math.sqrt(re*re + im*im);
				if (mod < seuil){
					FI.set_p_reel(i, j, 0);
					FI.set_p_imag(i, j, 0);
				} else{
					gardes++;
				}
			}
		}		
		return gardes;
	}

	public static CpxImg transformeCpxImg(CpxImg T) {
		double N = T.taille();
		double MaxPic = 0.0;
		double[][] log_L = new double[(int) N] [(int) N];
		double[][] F_log_L = new double[(int) N] [(int) N];
		CpxImg res = new CpxImg((int) N);

		for (int k = 0; k < N; k++) {
			for (int l = 0; l < N; l++) {
				double re = T.get_p_reel(k, l);
				double im = T.get_p_imag(k, l);
				double mod = Math.sqrt(re*re + im*im);
				double log = Math.log(1 + mod);
				log_L[k][l] = log;
				if (log > MaxPic) {
					MaxPic = log;
				}
			}
		}

		for (int k = 0; k < N; k++) {
			for (int l = 0; l < N; l++) {
				F_log_L[k][l] = (log_L[k][l] * 255)/MaxPic;
				res.set_p_reel(k, l, F_log_L[k][l]);
			}
		}
		return res;
	}

	static CpxImg annuler(CpxImg I, int k, int l){
		CpxImg TFD_I = FFT(I);
		TFD_I.set_p_reel(k, l, 0);
		TFD_I.set_p_imag(k, l, 0);
		CpxImg res = FFT_inverse(TFD_I);

		double erreurMax = 0.0;
			for (int i = 0; i < I.taille(); i++) {
				for (int j = 0; j < I.taille(); j++) {
					double diff = Math.abs(I.get_p_reel(i, j) - res.get_p_reel(i, j));
					if (diff > erreurMax) {
						erreurMax = diff;
					}
				}
			}
			System.out.println("Erreur maximale de reconstruction (" + k + ", " + l + ") : " + erreurMax);

		
		return res;
	}

	// Recherche le seuil par dichotomie
	public static double trouverSeuil(CpxImg I, int cible) {
		double bas = 0;
		double haut = 200000;
		double milieu = 0;

		for (int tour = 0; tour < 40; tour++) {
			milieu = (bas + haut) / 2;
			CpxImg F = FFT(I);
			int gardes = compression_seuil(F, milieu);
			if (gardes > cible) {
				bas = milieu;
			} else {
				haut = milieu;
			}
		}

		return haut;
	}


	
	public static void main(String[] args) {
		
		try {			
			//PLACEZ ICI VOS TESTS en 2D
			// lecture
			BytePixmap BP = new BytePixmap("AC_tp2_part2_donnees/tigre_512.pgm");
			// BytePixmap BP2 = new BytePixmap("AC_tp2_part2_donnees/tigre_512.pgm");
			// BytePixmap BP3 = new BytePixmap("AC_tp2_part2_donnees/tigre_512.pgm");
			CpxImg I = new CpxImg(BP);
			CpxImg TFD_I = FFT(I);

			// int comprSeuil = compression_seuil(TFD_I, 30);
			// System.out.println(comprSeuil);
			

			// for (int k=4; k <= 128; k*=2 ) {
			// 	CpxImg TFD_I = FFT(I);
			// 	compression(TFD_I, k);
			// 	CpxImg TFD_inverse = FFT_inverse(TFD_I);
			// 	BP = TFD_inverse.convert_to_BytePixmap();
			// 	String path = "resultats/ex5_" + k + "_centre.pgm";
			// 	BP.write(path);
				
			// }

			// int n = I.taille();
			// CpxImg annuleCentre = annuler(I, n/2, n/2);
			// CpxImg annulProcheCentre = annuler(I, n/2 + 1, n/2);
			// CpxImg annnulHF = annuler(I, n/2 + 100, n/2);
			// CpxImg TFD_I = FFT(I);
			// CpxImg TFD_inverse_I = FFT_inverse(TFD_I);
			// CpxImg TFD_trans_I = transformeCpxImg(TFD_I);
			
			// double erreurMax = 0.0;
			// for (int k = 0; k < I.taille(); k++) {
			// 	for (int l = 0; l < I.taille(); l++) {
			// 		double diff = Math.abs(I.get_p_reel(k, l) - TFD_inverse_I.get_p_reel(k, l));
			// 		if (diff > erreurMax) {
			// 			erreurMax = diff;
			// 		}
			// 	}
			// }
			// System.out.println("Erreur maximale de reconstruction : " + erreurMax);
			
			// Compare methode sur barbara et tigre
			BytePixmap BP_t = new BytePixmap("AC_tp2_part2_donnees/tigre_512.pgm");
			BytePixmap BP_b = new BytePixmap("AC_tp2_part2_donnees/barbara_512.pgm");
			CpxImg T = new CpxImg(BP_t);
			CpxImg B = new CpxImg(BP_b);
			CpxImg TFD_T = FFT(T);
			CpxImg TFD_B = FFT(B);

			compression(TFD_T, 81);
			compression(TFD_B, 81);
			CpxImg TFD_inverse_T = FFT_inverse(TFD_T);
			CpxImg TFD_inverse_B = FFT_inverse(TFD_B);
			BP_t = TFD_inverse_T.convert_to_BytePixmap();
			BP_b = TFD_inverse_B.convert_to_BytePixmap();
			BP_t.write("resultats/ex7_comprTigre.pgm");
			BP_b.write("resultats/ex7_comprBarbara.pgm");


			double seuil_T = trouverSeuil(T, 4 * 81 * 81);
			double seuil_B = trouverSeuil(B, 4 * 81 * 81);	
			CpxImg TFD_T2 = FFT(T);
			CpxImg TFD_B2 = FFT(B);		
			int gardes_t = compression_seuil(TFD_T2, seuil_T);
			int gardes_b = compression_seuil(TFD_B2, seuil_B);
			System.out.println(gardes_t);
			System.out.println(gardes_b);
			CpxImg TFD_inverse_T2 = FFT_inverse(TFD_T2);
			CpxImg TFD_inverse_B2 = FFT_inverse(TFD_B2);
			BP_t = TFD_inverse_T2.convert_to_BytePixmap();
			BP_b = TFD_inverse_B2.convert_to_BytePixmap();
			BP_t.write("resultats/ex7_comprSeuilTigre.pgm");
			BP_b.write("resultats/ex7_comprSeuilBarbara.pgm");



			// écriture

			// BP = TFD_I.convert_to_BytePixmap();
			// BP = TFD_trans_I.convert_to_BytePixmap();
			// BP = TFD_inverse_I.convert_to_BytePixmap();
			// BP = annuleCentre.convert_to_BytePixmap();
			// BP2 = annulProcheCentre.convert_to_BytePixmap();
			// BP3 = annnulHF.convert_to_BytePixmap();
			// BP.write("resultats/ex4_centre.pgm");
			// BP2.write("resultats/ex4_ProcheCentre.pgm");
			// BP3.write("resultats/ex4_HF.pgm");
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

}
