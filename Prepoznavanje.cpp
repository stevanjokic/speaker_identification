#include <iostream>
#include <limits>
#include <fstream>
#include <vector>
#include <string>
#include <sys/stat.h>
#include <math.h>
#include <windows.h>
#include <mmsystem.h>
#include <stdio.h>
#include <io.h>

using namespace std;

#define NSEC 10;
#define PI 3.14159265

HWAVEIN hWaveIn;
WAVEHDR WaveInHdr;
MMRESULT result;
WAVEFORMATEX pFormat;

void CheckMMIOError(DWORD code);
void SaveWaveFile(char* ime_govornika, int obuka);

int brojprozora;

vector<double> ucitavanjewavfajla(const char *imefajla)
{
	fstream podaci;
	podaci.open(imefajla, ios::in | ios::out | ios::binary);
	
	struct stat results;    // za dobijanje velicine fajla u bajtima
	stat(imefajla,&results);     // http://courses.cs.vt.edu/~cs2604/fall02/binio.html#stat
	//cout << "\nFajl " << imefajla << " sadrzi " << results.st_size << " bajtova" << "\n";

	/*for(int i=0; i<44; i++)
	{
		podaci.get();  // prolazak kroz 44 bajta zaglavlja
		                // 25,26,27 i 28. bajt predstavljaju ucestanost odabiranja (LitleEndian)
	}*/
	vector<int> referentni(4);referentni[0]=100;referentni[1]=97;referentni[2]=116;referentni[3]=97;vector<int> trenutni(4);
	trenutni[0]=podaci.get();trenutni[1]=podaci.get();trenutni[2]=podaci.get();trenutni[3]=podaci.get();int brojac=4;//brojac bajtova kroz koje se proslo
	while((trenutni[0]!=100) | (trenutni[1]!=97) | (trenutni[2]!=116) | (trenutni[3]!=97))
	{trenutni[0]=trenutni[1];trenutni[1]=trenutni[2];trenutni[2]=trenutni[3];trenutni[3]=podaci.get();brojac++;}
	//cout<<"data "<<brojac<<endl;
	for(int i=0;i<4;i++){podaci.get();} //4 bajta koji oznacavaju broj bajtova podataka, nakon cega slede podaci

	vector<double> odbirci ((results.st_size-brojac-4)/2); //vektor za cuvanje odbiraka signala
	for(int j=0; j<((results.st_size-brojac-4)/2); j++)
	{
	
		double vrednostodbirkauV;
	int b = podaci.get(); //vrednost prvog (nizeg) bajta                              
	//cout << b << "\n";                                                                 
	int a = podaci.get(); //vrednost drugog (viseg) bajta
	if (a>=0 && a<=127)   //vrednost odbirka je pozitivna ukoliko pocinje nulom, heksa vrednosti 0000h - 7FFF, odnosno vrednosti viseg bajta 00h - 7Fh
	{
		vrednostodbirkauV = ((a<<8) | b)/32767.0;  //drugi nizi bajt se pomera za 8 bita u levo na pravo mesto, mesto viseg bajta, zatim se primenjuje "ili" (sabiranje) po bitima za konacnu vrednost
		odbirci[j] = vrednostodbirkauV; //deli se sa 32767 da bi se vrednosti nalazile u opsegu -1 do 1
		//cout << j+1 << ". odbirak u [V] iznosi: " << vrednostodbirkauV << "\n"; 
	}

	else   //u suprotnom vrednost odbirka je negativna
	{
		int medjrez = ((a<<8) | b);  //predstava negativnog broja -B u komplementu dvojke se vrsi kao C2(B)=2^N-B, pri cemu je B pozitivan broj, medjrez je C2(B), N=16 (odbirci su zapisani pomocu 16 bita)
		int vredodbirka = medjrez - 65536;  // 2^N=2^16=65536, prava negativna vrednost odbirka na osnovu def. kompl. 2-ke
		//cout << j+1 << ". odbirak iznosi: " << vredodbirka << "\n";
		vrednostodbirkauV = vredodbirka/32768.0;
		odbirci[j] = vrednostodbirkauV;
		//cout << j+1 << ". odbirak u [V] iznosi: " << vrednostodbirkauV << "\n";
	}
	}
	return odbirci;
}

vector<double> uklanjanje_tisine(vector<double> odbirci)
{
	vector<double>odbirci1;
	
	int segment=600; //oko 27.21ms, duzina segmenta od oko 90ms = 2000 odbiraka
	int broj_segmenata=odbirci.size()/segment;
	double en_prvog_segmenta=0.0;
	for(int i=0;i<broj_segmenata;i++)
	{
		double en_segmenta=0.0;
		for(int j=0;j<segment;j++)
		{
			en_segmenta=en_segmenta+odbirci[i*segment+j]*odbirci[i*segment+j];
		}
		//cout << "Energija " << i << ". segmenta " << en_segmenta << endl;
		if(i==0) {en_prvog_segmenta=en_segmenta; /*cout << "Energija " << i+1 << ". segmenta " << en_segmenta << endl;*/}
		if(en_segmenta>5*en_prvog_segmenta)  //prvi segment u posmatranom signalu je uzet kao predstavnik tisine 
		{odbirci1.insert(odbirci1.end(),odbirci.begin()+i*segment,odbirci.begin()+(i+1)*segment);}
	}
	return odbirci1;
}

int osnovna_ucestanost(vector<double> odbirci)
{
	int visina_glasa;
	vector<double> abs_suma;
	double max=0.0;
	for(int i=80; i<250; i++)  //granice pomeraja za izracunavanje autokorelacija   (80-250)
	{
		double zbir=0.0;
		for(int j=0; j<(odbirci.size()-i); j++)
		{
			zbir=zbir+abs(odbirci[i+j]*odbirci[j]); //racunanje apsolutne vrednosti autokorelacije za i-tu vrednost pomeraja
		}
		if(zbir>=max)  //pomeraj pri kojem autokorelacija ima najvecu vrednost odgovara visini glasa
		{
			max=zbir;
			visina_glasa=i;
		}
	}
	return visina_glasa;
}

vector<double> kmfft(float N, vector<double> odbirci, vector<int> obrnutired,vector<float> prozor,vector<float> costab,vector<float> sintab)
{//Racuna kvadrate modula FFT za jedan prozor signala
	vector<double> privre(N);
	vector<double> privim(N);

	for (int i=0; i<N; i++)    //dobijanje prozora signala za koji se racuna odgovarajuci vektor obelezja
	   {
		   privre[i] = odbirci[obrnutired[i]] * prozor[obrnutired[i]]; //za realan signal vrsi se inicijalizacija samo re dela
		   privim[i] = 0;     // imaginarni deo se inicijalizuje na nulu posto je signal realan
	   }

	// DIT FFT radix2 u 2^n tacaka sa racunanjem modula koeficijenata u poslednjoj iteraciji
	
	vector<double> kmdfft (N);  // vektor za kvadrate modula DFT koeficijenata
	int brit; //broj iteracija FFT postupka
	brit = (log(N)/log(2.0));
	
	for (int k=1; k<=brit; k++)
	{
		double a=2.0;
		int korak;
		korak = pow(a,k);   // korak=2^k
		for (int i=0; i<=(N/korak-1); i++)
		{
			for (int j=0; j<=(korak/2-1); j++)
			{
				int arg=0;
				arg = j*N/korak;
				double pomocni1=0, pomocni2=0, pomocni3=0, pomocni4=0;

				// X(2^k*i+j) + / - (W u 2^k tacaka)^j*X(2^k*i+j+2^k-1) = X(2^k*i+j)  /  X(2^k*i+j+2^k-1)
				// (W u 2^k tacaka)^j = costab[j*N/2^k] - Im*sintab[j*N/2^k]

				pomocni1 = privre[korak*i+j];  //realni deo prvog sabirka
				pomocni2 = costab[arg]*privre[korak*i+j+korak/2] + sintab[arg]*privim[korak*i+j+korak/2]; //realni deo drugog sabirka
				pomocni3 = privim[korak*i+j];  //imaginarni deo prvog sabirka
				pomocni4 = costab[arg]*privim[korak*i+j+korak/2] - sintab[arg]*privre[korak*i+j+korak/2];  //imaginarni deo drugog sabirka
						
				privre[korak*i+j] = pomocni1 + pomocni2;
				privim[korak*i+j] = pomocni3 + pomocni4;

				privre[korak*i+j+korak/2] = pomocni1 - pomocni2;
				privim[korak*i+j+korak/2] = pomocni3 - pomocni4;

				if (k == brit)
				{
					kmdfft[korak*i+j] = privre[korak*i+j]*privre[korak*i+j] + privim[korak*i+j]*privim[korak*i+j];
					//cout << korak*i+j << ". " << "kmdfft: " << kmdfft[korak*i+j] << endl;   //prva polovina fft koeficjenata
					
					kmdfft[korak*i+j+korak/2] = privre[korak*i+j+korak/2]*privre[korak*i+j+korak/2] + privim[korak*i+j+korak/2]*privim[korak*i+j+korak/2];
					//cout << korak*i+j+korak/2 << ". " << "kmdfft: " << kmdfft[korak*i+j+korak/2] << endl;  // druga polovina fft koeficijenata, ispis se vrsi naizmenicno
				}
			}//j
		}//i
	}//k
	return kmdfft;
}

double spektralni_maksimum(float N, vector<double> kmdfft)
{
	double maksimum = -1.0;
	int ucestanost = 0;
	for (int i=0; i<=N/2; i++) //posto je DFT simetricna oko tacke N/2
	{
		if(kmdfft[i]>maksimum)
		{
			maksimum = kmdfft[i];
			ucestanost = i;
		}
	}
	cout << "\n Emax = " << maksimum << " na ucestanosti " << ucestanost << endl;
	return 0;
}

double sigm(double x)
{
	double s=1.0/(1.0+exp(-0.01*x));
	return s;
}

double sigm1(double x)
{
	double s1=1.0/(1.0+exp(-0.5*x));
	return s1;
}

vector<vector<double>> dodatna_obelezja(float N,int pomerajprozora,vector<double> odbirci,vector<int> obrnutired,vector<float> prozor,vector<float> costab,vector<float> sintab,int dodatno_obelezje1,int dodatno_obelezje2,int broj_dodatnih_obelezja)
{	
	vector<double> privre(N);
	vector<double> privim(N);
	vector<vector<double>> vektori_dodatnih_obelezja(brojprozora,vector<double> (3));
	vector<int> vektor_ucestanosti_maksimuma1(brojprozora);
	vector<int> vektor_ucestanosti_maksimuma2(brojprozora);
	vector<int> vektor_ucestanosti_maksimuma3(brojprozora);
	double najveci_maksimum_po_prozorima = -10000.0;
	
	for (int l=0; l<brojprozora; l++)  //pocetak petlje koja prolazi kroz sve prozore
	{
		int broj_nula=0;//promenljiva uvedena radi provere da li su u razmatranom prozoru svi odbirci nulte vrednosti
		//posto bi tada energija bilo kog skupa frekvencijskih komponenata tog prozora bila nula
		//ako bi se na tu nultu vrednost energije primenio logaritam dobio bi se nedefinisani broj -1.IND#
		//umesto toga na energije u takvim prozorima ne treba primenjivati logaritam nego takvim komponentama samo pridruziti nulu
		//odnosno dodatna obelezja za takve prozore signala imace nulte vrednosti
	   for (int i=0; i<N; i++)    //dobijanje prozora signala za koji se racuna odgovarajuci vektor obelezja
	   {
		   privre[i] = odbirci[l*pomerajprozora+obrnutired[i]] * prozor[obrnutired[i]]; //za realan signal vrsi se inicijalizacija samo re dela
		   privim[i] = 0;     // imaginarni deo se inicijalizuje na nulu posto je signal realan

		   if(privre[i]==0.0){broj_nula++;}
	   }

	   if(broj_nula==N) { for(int i=0;i<3;i++) {vektori_dodatnih_obelezja[l][i]=0.0;} }
	   else{


	// DIT FFT radix2 u 512 tacaka sa racunanjem modula koeficijenata u poslednjoj iteraciji
	
	   vector<double> kmdfft (N);  // vektor za kvadrate modula DFT koeficijenata
	   int brit; //broj iteracija FFT postupka
		brit = (log(N)/log(2.0));
	
		for (int k=1; k<=brit; k++)
		{
			double a=2.0;
			int korak;
			korak = pow(a,k);   // korak=2^k
			for (int i=0; i<=(N/korak-1); i++)
			{
				for (int j=0; j<=(korak/2-1); j++)
				{
					int arg=0;
					arg = j*N/korak;
					double pomocni1=0, pomocni2=0, pomocni3=0, pomocni4=0;

					// X(2^k*i+j) + / - (W u 2^k tacaka)^j*X(2^k*i+j+2^k-1) = X(2^k*i+j)  /  X(2^k*i+j+2^k-1)
					// (W u 2^k tacaka)^j = costab[j*N/2^k] - Im*sintab[j*N/2^k]

					pomocni1 = privre[korak*i+j];  //realni deo prvog sabirka
					pomocni2 = costab[arg]*privre[korak*i+j+korak/2] + sintab[arg]*privim[korak*i+j+korak/2]; //realni deo drugog sabirka
					pomocni3 = privim[korak*i+j];  //imaginarni deo prvog sabirka
					pomocni4 = costab[arg]*privim[korak*i+j+korak/2] - sintab[arg]*privre[korak*i+j+korak/2];  //imaginarni deo drugog sabirka
						
					privre[korak*i+j] = pomocni1 + pomocni2;
					privim[korak*i+j] = pomocni3 + pomocni4;

					privre[korak*i+j+korak/2] = pomocni1 - pomocni2;
					privim[korak*i+j+korak/2] = pomocni3 - pomocni4;

					if (k == brit)
					{
						kmdfft[korak*i+j] = privre[korak*i+j]*privre[korak*i+j] + privim[korak*i+j]*privim[korak*i+j];
						//cout << korak*i+j << ". " << "kmdfft: " << kmdfft[korak*i+j] << endl;   //prva polovina fft koeficjenata
					
						kmdfft[korak*i+j+korak/2] = privre[korak*i+j+korak/2]*privre[korak*i+j+korak/2] + privim[korak*i+j+korak/2]*privim[korak*i+j+korak/2];
						//cout << korak*i+j+korak/2 << ". " << "kmdfft: " << kmdfft[korak*i+j+korak/2] << endl;  // druga polovina fft koeficijenata, ispis se vrsi naizmenicno
					}// if
				}//j
			}//i
		}//k

		int norm_ucestanost1 = 0;
		//if(broj_dodatnih_obelezja==1 || broj_dodatnih_obelezja==2 || broj_dodatnih_obelezja==3)
		//{
		//odredjivanje dodatnog obelezja za razmatrani prozor signala
		double maksimum1 = -10000.0; double energija_u_okolini_max1 = 0.0;
		for(int i=0;i<N/2;i++)  //i=25
		{
			if(log(kmdfft[i])>maksimum1)
			{
				maksimum1 = log(kmdfft[i]); norm_ucestanost1 = i;
			}
		}
	
		//vektor_dodatnih_obelezja[l] = maksimum; // /((double)norm_ucestanost);
		vektor_ucestanosti_maksimuma1[l] = norm_ucestanost1;//ucestanost najvece komponente po prvom kriterijumu u razmatranom prozoru
		if(maksimum1>najveci_maksimum_po_prozorima) {najveci_maksimum_po_prozorima = maksimum1;}
		//vektor_dodatnih_obelezja[l]=log(kmdfft[norm_ucestanost1])+log(kmdfft[norm_ucestanost1-1])+log(kmdfft[norm_ucestanost1+1])+log(kmdfft[norm_ucestanost1-2])+log(kmdfft[norm_ucestanost1+2]); //uzimanje u obzir i okruzenja maksimuma
		if(norm_ucestanost1==0){vektor_ucestanosti_maksimuma1[l]=1;//pridruzuje joj se vrednost 1 da se ne bi delilo sa 0 pri racunanju tog dodatnog obelezja
			energija_u_okolini_max1 = maksimum1+log(kmdfft[norm_ucestanost1+1])*sigm(-1.0)+log(kmdfft[norm_ucestanost1+2])*sigm(-2.0);}
		if(norm_ucestanost1==1){energija_u_okolini_max1 = maksimum1+log(kmdfft[norm_ucestanost1-1])*sigm(-1.0)+log(kmdfft[norm_ucestanost1+1])*sigm(-1.0)+log(kmdfft[norm_ucestanost1+2])*sigm(-2.0);}
		if(norm_ucestanost1>1){
			energija_u_okolini_max1 = maksimum1+log(kmdfft[norm_ucestanost1-1])*sigm(-1.0)+log(kmdfft[norm_ucestanost1+1])*sigm(-1.0)+log(kmdfft[norm_ucestanost1-2])*sigm(-2.0)+log(kmdfft[norm_ucestanost1+2])*sigm(-2.0);}
		if(dodatno_obelezje1 == 1){
			vektori_dodatnih_obelezja[l][0]=energija_u_okolini_max1; }		
		//}

		int norm_ucestanost2 = 0;
		//if(broj_dodatnih_obelezja==2 || broj_dodatnih_obelezja==3)
		//{
		double maksimum2 = -10000; double energija_u_okolini_max2 = 0.0;
		for(int i=0;i<N/2;i++) //i=25,racunanje drugog po redu maksimuma u spektru
		{
			//if((i!=norm_ucestanost1)&&(i!=norm_ucestanost1-1)&&(i!=norm_ucestanost1+1)&&(i!=norm_ucestanost1-2)&&(i!=norm_ucestanost1+2))
			if((i<norm_ucestanost1-2)||(i>norm_ucestanost1+10)) //+10
			{
				if(log(kmdfft[i])>maksimum2)
				{
					maksimum2 = log(kmdfft[i]); norm_ucestanost2=i;
				}
			}
		}
		vektor_ucestanosti_maksimuma2[l] = norm_ucestanost2;//ucestanost najvece komponente po drugom kriterijumu u razmatranom prozoru
		if(norm_ucestanost2==0){energija_u_okolini_max2=maksimum2+log(kmdfft[norm_ucestanost2+1])*sigm(-1.0)+log(kmdfft[norm_ucestanost2+2])*sigm(-2.0);
		vektor_ucestanosti_maksimuma2[l]=1;}//da se ne bi delilo sa 0 pri racunanju ovog dodatnog obelezja
		if(norm_ucestanost2==1){energija_u_okolini_max2 = maksimum2+log(kmdfft[norm_ucestanost2-1])*sigm(-1.0)+log(kmdfft[norm_ucestanost2+1])*sigm(-1.0)+log(kmdfft[norm_ucestanost2+2])*sigm(-2.0);}
		if(norm_ucestanost2>1){
			energija_u_okolini_max2 = maksimum2+log(kmdfft[norm_ucestanost2-1])*sigm(-1.0)+log(kmdfft[norm_ucestanost2+1])*sigm(-1.0)+log(kmdfft[norm_ucestanost2-2])*sigm(-2.0)+log(kmdfft[norm_ucestanost2+2])*sigm(-2.0);}
		//vektor_dodatnih_obelezja[l]=log(kmdfft[norm_ucestanost2])+log(kmdfft[norm_ucestanost2-1])+log(kmdfft[norm_ucestanost2+1])+log(kmdfft[norm_ucestanost2-2])+log(kmdfft[norm_ucestanost2+2]);
		if(dodatno_obelezje1==0 && dodatno_obelezje2==1) {vektori_dodatnih_obelezja[l][0]=energija_u_okolini_max2;}
		
		if(dodatno_obelezje1==1 && dodatno_obelezje2==1) {vektori_dodatnih_obelezja[l][1]=energija_u_okolini_max2;}
		  //vektor_dodatnih_obelezja[l]=(energ_u_okolini_max-energ_u_okolini_max2)/((double)(norm_ucestanost1-norm_ucestanost));
		//}izbacuju se ovi uslovi da bi se moglo eksperimentisati i samo sa postojanjem drugog dodatnog obelezja kada se ne koristi prvo dodatno obelezje

		int norm_ucestanost3 = 2;
		if(broj_dodatnih_obelezja==3)
		{
		double maksimum3 = -10000; double energija_u_okolini_max3 = 0.0;
		for(int i=0;i<N/2;i++) //i=25,racunanje treceg po redu maksimuma u spektru
		{
			if(((i<norm_ucestanost2-2)||(i>norm_ucestanost2+10)) && i!=norm_ucestanost1-3 && i!=norm_ucestanost1-2 && i!=norm_ucestanost1-1 &&
				i!=norm_ucestanost1 && i!=norm_ucestanost1+1 && i!=norm_ucestanost1+2 && i!=norm_ucestanost1+3 && i!=norm_ucestanost1+4 &&
				i!=norm_ucestanost1+5 && i!=norm_ucestanost1+6 && i!=norm_ucestanost1+7 && i!=norm_ucestanost1+8 && i!=norm_ucestanost1+9 &&
				i!=norm_ucestanost1+10) //+10
			{
				if(log(kmdfft[i])>maksimum3)
				{
					maksimum3 = log(kmdfft[i]); norm_ucestanost3=i;
				}
			}
		}
		vektor_ucestanosti_maksimuma3[l] = norm_ucestanost3;//ucestanost najvece komponente po trecem kriterijumu u razmatranom prozoru
		if(norm_ucestanost3==0){energija_u_okolini_max3=maksimum3+log(kmdfft[norm_ucestanost3+1])*sigm(-1.0)+log(kmdfft[norm_ucestanost3+2])*sigm(-2.0);//*sigm(-1.0)//*sigm(-2.0)
		vektor_ucestanosti_maksimuma3[l]=1;}//da se ne bi delilo sa 0 pri racunanju ovog dodatnog obelezja
		if(norm_ucestanost3==1){energija_u_okolini_max3 = maksimum3+log(kmdfft[norm_ucestanost3-1])*sigm(-1.0)+log(kmdfft[norm_ucestanost3+1])*sigm(-1.0)+log(kmdfft[norm_ucestanost3+2])*sigm(-2.0);}
		if(norm_ucestanost3>1){
			energija_u_okolini_max3 = maksimum3+log(kmdfft[norm_ucestanost3-1])*sigm(-1.0)+log(kmdfft[norm_ucestanost3+1])*sigm(-1.0)+log(kmdfft[norm_ucestanost3-2])*sigm(-2.0)+log(kmdfft[norm_ucestanost3+2])*sigm(-2.0);}
		vektori_dodatnih_obelezja[l][2]=energija_u_okolini_max3;
		}
		}//kraj od else
	}//l

	//ponderisanje dodatnih obelezja shodno njihovim vrednostima u odnosu na najveci maksimum po prozorima
	for(int i=0;i<brojprozora;i++)
	{
		if(vektori_dodatnih_obelezja[i][0]==0.0 && vektori_dodatnih_obelezja[i][1]==0.0 && vektori_dodatnih_obelezja[i][2]==0.0)
		{vektori_dodatnih_obelezja[i][0]=0.0;vektori_dodatnih_obelezja[i][1]=0.0;vektori_dodatnih_obelezja[i][2]=0.0;}
		else{
		//if(broj_dodatnih_obelezja==1 || broj_dodatnih_obelezja==2 || broj_dodatnih_obelezja==3)
		if(dodatno_obelezje1==1)
		{vektori_dodatnih_obelezja[i][0] = (vektori_dodatnih_obelezja[i][0]*vektori_dodatnih_obelezja[i][0]/najveci_maksimum_po_prozorima)/((double)vektor_ucestanosti_maksimuma1[i]);
		//int k1=vektor_ucestanosti_maksimuma1[i]; //int k1d=k1-2, k1g=k1+2; double f1d=44100.0*((double)k1d), f1g=44100.0*((double)k1g);
		//double f1=44100.0*((double)k1), f1_mel=2595.0*log(1.0+f1/700.0); //f1d_mel=2595.0*log(1.0+f1d/700.0), f1g_mel=2595.0*log(1.0+f1g/700.0);
		//double df1_mel=f1g_mel-f1d_mel;
		//vektori_dodatnih_obelezja[i][0] = (vektori_dodatnih_obelezja[i][0]*vektori_dodatnih_obelezja[i][0]/najveci_maksimum_po_prozorima)/f1_mel;
		}
		//if(broj_dodatnih_obelezja==2 || broj_dodatnih_obelezja==3)
		if(dodatno_obelezje1==0 && dodatno_obelezje2==1)
		{vektori_dodatnih_obelezja[i][0] = (vektori_dodatnih_obelezja[i][0]*vektori_dodatnih_obelezja[i][0]/najveci_maksimum_po_prozorima)/((double)vektor_ucestanosti_maksimuma2[i]);}
		if(dodatno_obelezje1==1 && dodatno_obelezje2==1)
		{vektori_dodatnih_obelezja[i][1] = (vektori_dodatnih_obelezja[i][1]*vektori_dodatnih_obelezja[i][1]/najveci_maksimum_po_prozorima)/((double)vektor_ucestanosti_maksimuma2[i]);
		//int k2=vektor_ucestanosti_maksimuma2[i];
		//double f2=44100.0*((double)k2), f2_mel=2595.0*log(1.0+f2/700.0);
		//vektori_dodatnih_obelezja[i][1] = (vektori_dodatnih_obelezja[i][1]*vektori_dodatnih_obelezja[i][1]/najveci_maksimum_po_prozorima)/f2_mel;
		}
		if(broj_dodatnih_obelezja==3)
		{vektori_dodatnih_obelezja[i][2] = (vektori_dodatnih_obelezja[i][2]*vektori_dodatnih_obelezja[i][2]/najveci_maksimum_po_prozorima)/((double)vektor_ucestanosti_maksimuma3[i]);}
		}//kraj od else
	}
	return vektori_dodatnih_obelezja;
}

vector<vector<double>> mfcc(float N,int pomerajprozora,vector<double> odbirci,char *oblik_kr_opsega,int brkropsega,int brojMFCCs,vector<int> obrnutired,vector<float> prozor,vector<float> costab,vector<float> sintab,vector<vector<float>> costabMFCC)
{	
	vector<double> privre(N);
	vector<double> privim(N);
	vector<vector<double>> vektoriobelezja(brojprozora, vector<double> (brojMFCCs));
	
	for (int l=0; l<brojprozora; l++)  //pocetak petlje koja prolazi kroz sve prozore
	{
	   for (int i=0; i<N; i++)    //dobijanje prozora signala za koji se racuna odgovarajuci vektor obelezja
	   {
		   privre[i] = odbirci[l*pomerajprozora+obrnutired[i]] * prozor[obrnutired[i]]; //za realan signal vrsi se inicijalizacija samo re dela
		   privim[i] = 0;     // imaginarni deo se inicijalizuje na nulu posto je signal realan
	   }


	// DIT FFT radix2 u 512 tacaka sa racunanjem modula koeficijenata u poslednjoj iteraciji
	
	vector<double> kmdfft (N);  // vektor za kvadrate modula DFT koeficijenata
	int brit; //broj iteracija FFT postupka
	brit = (log(N)/log(2.0));
	
	for (int k=1; k<=brit; k++)
	{
		double a=2.0;
		int korak;
		korak = pow(a,k);   // korak=2^k
		for (int i=0; i<=(N/korak-1); i++)
		{
			for (int j=0; j<=(korak/2-1); j++)
			{
				int arg=0;
				arg = j*N/korak;
				double pomocni1=0, pomocni2=0, pomocni3=0, pomocni4=0;

				// X(2^k*i+j) + / - (W u 2^k tacaka)^j*X(2^k*i+j+2^k-1) = X(2^k*i+j)  /  X(2^k*i+j+2^k-1)
				// (W u 2^k tacaka)^j = costab[j*N/2^k] - Im*sintab[j*N/2^k]

				pomocni1 = privre[korak*i+j];  //realni deo prvog sabirka
				pomocni2 = costab[arg]*privre[korak*i+j+korak/2] + sintab[arg]*privim[korak*i+j+korak/2]; //realni deo drugog sabirka
				pomocni3 = privim[korak*i+j];  //imaginarni deo prvog sabirka
				pomocni4 = costab[arg]*privim[korak*i+j+korak/2] - sintab[arg]*privre[korak*i+j+korak/2];  //imaginarni deo drugog sabirka
						
				privre[korak*i+j] = pomocni1 + pomocni2;
				privim[korak*i+j] = pomocni3 + pomocni4;

				privre[korak*i+j+korak/2] = pomocni1 - pomocni2;
				privim[korak*i+j+korak/2] = pomocni3 - pomocni4;

				if (k == brit)
				{
					kmdfft[korak*i+j] = privre[korak*i+j]*privre[korak*i+j] + privim[korak*i+j]*privim[korak*i+j];
					//cout << korak*i+j << ". " << "kmdfft: " << kmdfft[korak*i+j] << endl;   //prva polovina fft koeficjenata
					
					kmdfft[korak*i+j+korak/2] = privre[korak*i+j+korak/2]*privre[korak*i+j+korak/2] + privim[korak*i+j+korak/2]*privim[korak*i+j+korak/2];
					//cout << korak*i+j+korak/2 << ". " << "kmdfft: " << kmdfft[korak*i+j+korak/2] << endl;  // druga polovina fft koeficijenata, ispis se vrsi naizmenicno
				}
			}
		}
	}
	
	//racunanje log energija u odgovarajucim kriticnim opsezima, korisceni su pravougaoni filtri, 20 opsega
	vector<double> logenergijaopsega (brkropsega);//vektor za cuvanje log energija posmatranog l-tog prozora signala
	int i; //redni broj filtra
	int k,k1,k2;
	float s,s1;
	float zbir;
	float a,b,c=10.0, fs=22050.0;
	for (i=0; i<brkropsega; i++)
	{
		zbir = 0.0;
		s = i*150.0/2595.0;
		a = pow(c,s);
		k1 = N * (a - 1.0) * 700.0/fs;  //donja granica i-tog opsega
		s1 = (i+2)*150.0/2595.0;
		b = pow(c,s1);
		k2 = N * (b - 1.0) * 700.0/fs; //gornja granica i-tog opsega
		
		double kc=(k1+k2)/2.0;  //srediste kriticnog opsega u Hz skali
		//float sc=(i+1)*100.0/2595.0;  //eksperimentalno srediste
		//float sc=(s+s1)/2;  //srediste kriticnog opsega u mel skali
		//float ac=pow(c,sc); int kc=N*(ac-1.0)*700.0/22050;  //centralna ucestanost kriticnog opsega

		/*int kc=k1;  //maksimum kriticnog opsega tamo gde je najveca komponenta
		for(k=k1; k<=k2; k++)
		{
			if(kmdfft[k]>kmdfft[kc]) {kc=k;}
		}*/
		
		string tip; tip=string(oblik_kr_opsega); //potrebno radi konverzije *char u string, da bi se poredili nizovi znakova potrebno ih je definisati kao string tipove 
		for (k=k1; k<=k2; k++)
		{
			if(tip=="p") {zbir += (kmdfft[k]);}   //pravougaoni kriticni opseg
			
			if(tip=="t")
			{if(k>=k1 && k<=(k2+k1)/2.0) {zbir=zbir+2.0/(k2-k1)*(k-k1)*kmdfft[k];}  //trougaoni kriticni opseg, prva strmina
			if(k>(k2+k1)/2.0 && k<=k2) {zbir=zbir+2.0/(k1-k2)*(k-k2)*kmdfft[k];}}   //druga strmina
			
			//if(k>=k1 && k<=kc) {zbir=zbir+((float)exp((float)((k-k1)*1.0))-1.0)*kmdfft[k];} //eksponencijalni kriticni opseg sa prilagodljivim sredistem (gornji deo eksp funkcije), prva strmina
			//if(k>kc && k<=k2) {zbir=zbir+((float)exp((float)((k2-k)*1.0))-1.0)*kmdfft[k];}
			
			if(tip=="e_g-1")
			{if(k>=k1 && k<=kc) {zbir=zbir+((float)exp((float)((k-k1)*1.0))-1.0)*kmdfft[k];} //eksponencijalni kriticni opseg, prva strmina (gornji deo eksp funkcije) (spusten za 1)
			if(k>kc && k<=k2) {zbir=zbir+((float)exp((float)((k2-k)*1.0))-1.0)*kmdfft[k];}}

			if(tip=="e_g")
			{if(k>=k1 && k<=kc) {zbir=zbir+((float)exp((float)((k-k1)*1.0)))*kmdfft[k];} //eksponencijalni kriticni opseg, prva strmina (gornji deo eksp funkcije)
			if(k>kc && k<=k2) {zbir=zbir+((float)exp((float)((k2-k)*1.0)))*kmdfft[k];}}

			//float argument1=k-k1;    // 2*PI*(k-k1)/(kc-k1);
			//float argument1=kc-k;
			//float argument2=k2-k;     //2*PI*(k2-k)/(k2-kc);
			//float argument2=k-kc;
			//if(k>=k1 && k<=kc) {zbir=zbir+((float)exp((float)((k-k1)*1.0)))*kmdfft[k]*abs((float)cos(argument1));} //eksponencijalni-cos kriticni opseg, prva strmina (gornji deo eksp funkcije)
			//if(k>kc && k<=k2) {zbir=zbir+((float)exp((float)((k2-k)*1.0)))*kmdfft[k]*abs((float)cos(argument2));}

			if(tip=="e_d")
			{if(k>=k1 && k<=kc) {zbir=zbir+((float)exp((float)((k-kc)*2.0)))*kmdfft[k];} //eksponencijalni kriticni opseg, prva strmina (donji deo eksp funkcije)
			if(k>kc && k<=k2) {zbir=zbir+((float)exp((float)((kc-k)*2.0)))*kmdfft[k];}}

			if(tip=="sigm")
			{if(k>=k1 && k<kc) {zbir=zbir+(sigm1((double)(k)-(double)(kc)))*kmdfft[k];}
			if(k==kc) {zbir=zbir+kmdfft[k];}
			if(k>kc && k<=k2) {zbir=zbir+(sigm1((double)(kc)-(double)(k)))*kmdfft[k];}}

			//float osnova=4.0;
			//if(k>=k1 && k<=(k2+k1)/2.0) {zbir=zbir+pow(osnova,(float)2.0*(k-k1))*kmdfft[k];}
			//if(k>(k2+k1)/2.0 && k<=k2) {zbir=zbir+pow(osnova,(float)2.0*(k2-k))*kmdfft[k];}
			
			  //if(kmdfft[k]>max_kmdfft) {max_kmdfft=kmdfft[k];} 
			
			/*float br_podopsega=40.0;
			for(int n=0;n<br_podopsega;n++) 
			{ 
				if ((k>=k1+n*(k2-k1)/br_podopsega) && (k<=k1+(k2-k1)/br_podopsega/2.0+n*(k2-k1)/br_podopsega))  //prva polovina podopsega
				{zbir=zbir+((float)exp((float)((k-k1-n*(k2-k1)/br_podopsega)*2.0)))*kmdfft[k];}
				if ((k>k1+(k2-k1)/br_podopsega/2.0+n*(k2-k1)/br_podopsega) && (k<k1+(n+1)*(k2-k1)/br_podopsega))
				{zbir=zbir+((float)exp((float)((k1+(n+1)*(k2-k1)/br_podopsega-k)*2.0)))*kmdfft[k];}
			}*/
		}
		//zbir=max_kmdfft; //uzimanje samo fft koeficijenta koji ima najvecu energiju u posmatranom kriticnom opsegu
		if(zbir>0) {logenergijaopsega[i] = (log(2.0*zbir)/log(10.0));}
		if(zbir=0) {logenergijaopsega[i]=0;} //moze se desiti dugacak niz nula u snimku, nizovi ramova da budu nule,
		//verovatno iz nekih razloga da tada nije nista snimljeno, rezultat bi bio log(0) (-1.#IND)
	}
	
	double min_log_energ=logenergijaopsega[0]; //pretpostavka da najmanja log_energ predstavlja sum, koji je onda bolje izdvojiti radi smanjenja njegovog uticaja na MFCCs
	vector<double> vaznost_kriticnih_opsega(brkropsega); for(int i=0;i<brkropsega;i++) {vaznost_kriticnih_opsega[i]=1.0;} //postavljanje pocetnih vrednosti koeficijenata vaznosti
	for(int i=0;i<brkropsega;i++) { if(logenergijaopsega[i]<min_log_energ) {min_log_energ=logenergijaopsega[i];}} //odredjivanje najmanje logenerg u posmatranim kriticnim opsezima
	//for(int i=0;i<brkropsega;i++) { if(logenergijaopsega[i]<=min_log_energ) {logenergijaopsega[i]=0.0;}} //ova popravka nad CHAINS bazom dosta pogorsava tacnost
	//for(int i=0;i<brkropsega;i++) {logenergijaopsega[i]=logenergijaopsega[i]-min_log_energ;} //ova popravka nad CHAINS bazom kada se koristi po jedan fajl za obuku i za test malo pogorsa tacnost

	//odredjivanje MFCCs
	
		for (int j=0; j<brojMFCCs; j++) 
		{
			for (int i=0; i<brkropsega; i++)  
			{
				vektoriobelezja[l] [j] += logenergijaopsega[i]*costabMFCC[j][i];
			}
		}
	}
	return vektoriobelezja;
}

vector<vector<double>> obelezja(float N,int pomerajprozora,vector<double> odbirci,vector<int> obrnutired,vector<float> prozor,vector<float> costab,vector<float> sintab,int broj_obelezja)
{	
	vector<double> privre(N);
	vector<double> privim(N);
	vector<vector<double>> vekt_obelezja(brojprozora,vector<double> (broj_obelezja));
	vector<int> vektor_ucestanosti_maksimuma1(brojprozora);
	vector<int> vektor_ucestanosti_maksimuma2(brojprozora);
	vector<int> vektor_ucestanosti_maksimuma3(brojprozora);
	vector<int> vektor_ucestanosti_maksimuma4(brojprozora);
	vector<int> vektor_ucestanosti_maksimuma5(brojprozora);
	vector<int> vektor_ucestanosti_maksimuma6(brojprozora);
	vector<int> vektor_ucestanosti_maksimuma7(brojprozora);
	vector<int> vektor_ucestanosti_maksimuma8(brojprozora);
	vector<int> vektor_ucestanosti_maksimuma9(brojprozora);
	vector<int> vektor_ucestanosti_maksimuma10(brojprozora);
	vector<int> vektor_ucestanosti_maksimuma11(brojprozora);
	vector<int> vektor_ucestanosti_maksimuma12(brojprozora);
	vector<vector<int>> vektor_ucestanosti_maksimuma(brojprozora, vector<int> (broj_obelezja));
	
	vector<vector<double>> log_energije(brojprozora,vector<double> (broj_obelezja));
	
	double najveci_maksimum_po_prozorima = -10000.0;
	
	for (int l=0; l<brojprozora; l++)  //pocetak petlje koja prolazi kroz sve prozore
	{
	   for (int i=0; i<N; i++)    //dobijanje prozora signala za koji se racuna odgovarajuci vektor obelezja
	   {
		   privre[i] = odbirci[l*pomerajprozora+obrnutired[i]] * prozor[obrnutired[i]]; //za realan signal vrsi se inicijalizacija samo re dela
		   privim[i] = 0;     // imaginarni deo se inicijalizuje na nulu posto je signal realan
	   }


	// DIT FFT radix2 u 512 tacaka sa racunanjem modula koeficijenata u poslednjoj iteraciji
	
	   vector<double> kmdfft (N);  // vektor za kvadrate modula DFT koeficijenata
	   vector<double> sigmkmdfft (N);  //primenjena sigmoid funkcija na kmdfft
	   int brit; //broj iteracija FFT postupka
		brit = (log(N)/log(2.0));
	
		for (int k=1; k<=brit; k++)
		{
			double a=2.0;
			int korak;
			korak = pow(a,k);   // korak=2^k
			for (int i=0; i<=(N/korak-1); i++)
			{
				for (int j=0; j<=(korak/2-1); j++)
				{
					int arg=0;
					arg = j*N/korak;
					double pomocni1=0, pomocni2=0, pomocni3=0, pomocni4=0;

					// X(2^k*i+j) + / - (W u 2^k tacaka)^j*X(2^k*i+j+2^k-1) = X(2^k*i+j)  /  X(2^k*i+j+2^k-1)
					// (W u 2^k tacaka)^j = costab[j*N/2^k] - Im*sintab[j*N/2^k]

					pomocni1 = privre[korak*i+j];  //realni deo prvog sabirka
					pomocni2 = costab[arg]*privre[korak*i+j+korak/2] + sintab[arg]*privim[korak*i+j+korak/2]; //realni deo drugog sabirka
					pomocni3 = privim[korak*i+j];  //imaginarni deo prvog sabirka
					pomocni4 = costab[arg]*privim[korak*i+j+korak/2] - sintab[arg]*privre[korak*i+j+korak/2];  //imaginarni deo drugog sabirka
						
					privre[korak*i+j] = pomocni1 + pomocni2;
					privim[korak*i+j] = pomocni3 + pomocni4;

					privre[korak*i+j+korak/2] = pomocni1 - pomocni2;
					privim[korak*i+j+korak/2] = pomocni3 - pomocni4;

					if (k == brit)
					{
						kmdfft[korak*i+j] = privre[korak*i+j]*privre[korak*i+j] + privim[korak*i+j]*privim[korak*i+j];
						//cout << korak*i+j << ". " << "kmdfft: " << kmdfft[korak*i+j] << endl;   //prva polovina fft koeficijenata
						sigmkmdfft[korak*i+j] = 1.0/(1.0+exp(0.0-kmdfft[korak*i+j]));
					
						kmdfft[korak*i+j+korak/2] = privre[korak*i+j+korak/2]*privre[korak*i+j+korak/2] + privim[korak*i+j+korak/2]*privim[korak*i+j+korak/2];
						//cout << korak*i+j+korak/2 << ". " << "kmdfft: " << kmdfft[korak*i+j+korak/2] << endl;  // druga polovina fft koeficijenata, ispis se vrsi naizmenicno
						sigmkmdfft[korak*i+j+korak/2] = 1.0/(1.0+exp(0.0-kmdfft[korak*i+j+korak/2]));
					}// if
				}//j
			}//i
		}//k

		int norm_ucestanost1 = 0;
		if(broj_obelezja==1 || broj_obelezja==2 || broj_obelezja==3 || broj_obelezja==4 || broj_obelezja==5 || broj_obelezja==6 || broj_obelezja==7 || broj_obelezja==8 || broj_obelezja==9 || broj_obelezja==10 || broj_obelezja==11 || broj_obelezja==12)
		{
		//odredjivanje obelezja za razmatrani prozor signala
		double maksimum1 = -10000.0; double energija_u_okolini_max1 = 0.0;
		for(int i=1;i<N/2;i++)
		{
			/*if(log(kmdfft[i])>maksimum1)
			{
				maksimum1 = log(kmdfft[i]); norm_ucestanost1 = i;
			}*/
			if(sigmkmdfft[i]>maksimum1)
			{
				maksimum1 = sigmkmdfft[i]; norm_ucestanost1 = i;
			}
		}
	
		//vektor_dodatnih_obelezja[l] = maksimum; // /((double)norm_ucestanost);
		vektor_ucestanosti_maksimuma1[l] = norm_ucestanost1;//ucestanost najvece komponente po prvom kriterijumu u razmatranom prozoru
		vektor_ucestanosti_maksimuma[l][0]=norm_ucestanost1;
		if(maksimum1>najveci_maksimum_po_prozorima) {najveci_maksimum_po_prozorima = maksimum1;}
		//vektor_dodatnih_obelezja[l]=log(kmdfft[norm_ucestanost1])+log(kmdfft[norm_ucestanost1-1])+log(kmdfft[norm_ucestanost1+1])+log(kmdfft[norm_ucestanost1-2])+log(kmdfft[norm_ucestanost1+2]); //uzimanje u obzir i okruzenja maksimuma
		if(norm_ucestanost1>3)
		//{energija_u_okolini_max1=maksimum1+log(kmdfft[norm_ucestanost1-1])+log(kmdfft[norm_ucestanost1+1])+log(kmdfft[norm_ucestanost1-2])+log(kmdfft[norm_ucestanost1+2])+log(kmdfft[norm_ucestanost1-3])+log(kmdfft[norm_ucestanost1+3])+log(kmdfft[norm_ucestanost1-4])+log(kmdfft[norm_ucestanost1+4]);}
		{energija_u_okolini_max1=maksimum1+sigmkmdfft[norm_ucestanost1-1]+sigmkmdfft[norm_ucestanost1+1]+sigmkmdfft[norm_ucestanost1-2]+sigmkmdfft[norm_ucestanost1+2]+sigmkmdfft[norm_ucestanost1-3]+sigmkmdfft[norm_ucestanost1+3]+sigmkmdfft[norm_ucestanost1-4]+sigmkmdfft[norm_ucestanost1+4];}
		if(norm_ucestanost1>2)
		//{energija_u_okolini_max1=maksimum1+log(kmdfft[norm_ucestanost1-1])+log(kmdfft[norm_ucestanost1+1])+log(kmdfft[norm_ucestanost1-2])+log(kmdfft[norm_ucestanost1+2])+log(kmdfft[norm_ucestanost1-3])+log(kmdfft[norm_ucestanost1+3]);}
		{energija_u_okolini_max1=maksimum1+sigmkmdfft[norm_ucestanost1-1]+sigmkmdfft[norm_ucestanost1+1]+sigmkmdfft[norm_ucestanost1-2]+sigmkmdfft[norm_ucestanost1+2]+sigmkmdfft[norm_ucestanost1-3]+sigmkmdfft[norm_ucestanost1+3];}
		if(norm_ucestanost1>1)
		//{energija_u_okolini_max1 = maksimum1+log(kmdfft[norm_ucestanost1-1])+log(kmdfft[norm_ucestanost1+1])+log(kmdfft[norm_ucestanost1-2])+log(kmdfft[norm_ucestanost1+2])+log(kmdfft[norm_ucestanost1+3])+log(kmdfft[norm_ucestanost1+4])+log(kmdfft[norm_ucestanost1+5]);}
		{energija_u_okolini_max1 = maksimum1+sigmkmdfft[norm_ucestanost1-1]+sigmkmdfft[norm_ucestanost1+1]+sigmkmdfft[norm_ucestanost1-2]+sigmkmdfft[norm_ucestanost1+2]+sigmkmdfft[norm_ucestanost1+3]+sigmkmdfft[norm_ucestanost1+4]+sigmkmdfft[norm_ucestanost1+5];}
		if(norm_ucestanost1==1)
		//{energija_u_okolini_max1 = maksimum1+log(kmdfft[norm_ucestanost1-1])+log(kmdfft[norm_ucestanost1+1])+log(kmdfft[norm_ucestanost1+2])+log(kmdfft[norm_ucestanost1+3])+log(kmdfft[norm_ucestanost1+4])+log(kmdfft[norm_ucestanost1+5]);}
		{energija_u_okolini_max1 = maksimum1+sigmkmdfft[norm_ucestanost1-1]+sigmkmdfft[norm_ucestanost1+1]+sigmkmdfft[norm_ucestanost1+2]+sigmkmdfft[norm_ucestanost1+3]+sigmkmdfft[norm_ucestanost1+4]+sigmkmdfft[norm_ucestanost1+5];}
		if(norm_ucestanost1==0)
		//{energija_u_okolini_max1 = maksimum1+log(kmdfft[norm_ucestanost1+1])+log(kmdfft[norm_ucestanost1+2])+log(kmdfft[norm_ucestanost1+3])+log(kmdfft[norm_ucestanost1+4])+log(kmdfft[norm_ucestanost1+5]);}
		{energija_u_okolini_max1 = maksimum1+sigmkmdfft[norm_ucestanost1+1]+sigmkmdfft[norm_ucestanost1+2]+sigmkmdfft[norm_ucestanost1+3]+sigmkmdfft[norm_ucestanost1+4]+sigmkmdfft[norm_ucestanost1+5];}
		vekt_obelezja[l][0]=energija_u_okolini_max1; log_energije[l][0]=energija_u_okolini_max1;		
		}

		int norm_ucestanost2 = 0;
		if(broj_obelezja==2 || broj_obelezja==3 || broj_obelezja==4 || broj_obelezja==5 || broj_obelezja==6 || broj_obelezja==7 || broj_obelezja==8 || broj_obelezja==9 || broj_obelezja==10 || broj_obelezja==11 || broj_obelezja==12)
		{
		double maksimum2 = -10000; double energija_u_okolini_max2 = 0.0;
		for(int i=1;i<N/2;i++) //racunanje drugog po redu maksimuma u spektru
		{
			//if((i!=norm_ucestanost1)&&(i!=norm_ucestanost1-1)&&(i!=norm_ucestanost1+1)&&(i!=norm_ucestanost1-2)&&(i!=norm_ucestanost1+2))
			if((i<norm_ucestanost1-5) || (i>norm_ucestanost1+5))
			{
				if(log(kmdfft[i])>maksimum2)
				{
					maksimum2 = log(kmdfft[i]); norm_ucestanost2=i;
				}
			}
		}
		vektor_ucestanosti_maksimuma2[l] = norm_ucestanost2;//ucestanost najvece komponente po drugom kriterijumu u razmatranom prozoru
		vektor_ucestanosti_maksimuma[l][1]=norm_ucestanost2;
		if(norm_ucestanost2>2)
		{energija_u_okolini_max2=maksimum2+log(kmdfft[norm_ucestanost2-1])+log(kmdfft[norm_ucestanost2+1])+log(kmdfft[norm_ucestanost2-2])+log(kmdfft[norm_ucestanost2+2])+log(kmdfft[norm_ucestanost2-3])+log(kmdfft[norm_ucestanost2+3]);}
		if(norm_ucestanost2>1)
		{energija_u_okolini_max2 = maksimum2+log(kmdfft[norm_ucestanost2-1])+log(kmdfft[norm_ucestanost2+1])+log(kmdfft[norm_ucestanost2-2])+log(kmdfft[norm_ucestanost2+2]);}
		if(norm_ucestanost2==1)
		{energija_u_okolini_max2 = maksimum2+log(kmdfft[norm_ucestanost2-1])+log(kmdfft[norm_ucestanost2+1])+log(kmdfft[norm_ucestanost2+2]);}
		if(norm_ucestanost2==0)
		{energija_u_okolini_max2 = maksimum2+log(kmdfft[norm_ucestanost2+1])+log(kmdfft[norm_ucestanost2+2]);}
		//vektor_dodatnih_obelezja[l]=log(kmdfft[norm_ucestanost2])+log(kmdfft[norm_ucestanost2-1])+log(kmdfft[norm_ucestanost2+1])+log(kmdfft[norm_ucestanost2-2])+log(kmdfft[norm_ucestanost2+2]);
		vekt_obelezja[l][1]=energija_u_okolini_max2; log_energije[l][1]=energija_u_okolini_max2;
		  //vektor_dodatnih_obelezja[l]=(energ_u_okolini_max-energ_u_okolini_max2)/((double)(norm_ucestanost1-norm_ucestanost));
		}

		int norm_ucestanost3 = 2;
		if(broj_obelezja==3 || broj_obelezja==4 || broj_obelezja==5 || broj_obelezja==6 || broj_obelezja==7 || broj_obelezja==8 || broj_obelezja==9 || broj_obelezja==10 || broj_obelezja==11 || broj_obelezja==12)
		{
		double maksimum3 = -10000; double energija_u_okolini_max3 = 0.0;
		for(int i=1;i<N/2;i++) //racunanje treceg po redu maksimuma u spektru
		{
			if( ((i<norm_ucestanost1-5) || (i>norm_ucestanost1+5)) && ((i<norm_ucestanost2-5) || (i>norm_ucestanost2+5)) )
			{
				if(log(kmdfft[i])>maksimum3)
				{
					maksimum3 = log(kmdfft[i]); norm_ucestanost3=i;
				}
			}
		}
		vektor_ucestanosti_maksimuma3[l] = norm_ucestanost3;//ucestanost najvece komponente po trecem kriterijumu u razmatranom prozoru
		vektor_ucestanosti_maksimuma[l][2]=norm_ucestanost3;
		if(norm_ucestanost3>3)
		{energija_u_okolini_max3 = maksimum3+log(kmdfft[norm_ucestanost3-1])+log(kmdfft[norm_ucestanost3+1])+log(kmdfft[norm_ucestanost3-2])+log(kmdfft[norm_ucestanost3+2])+log(kmdfft[norm_ucestanost3-3])+log(kmdfft[norm_ucestanost3+3])+log(kmdfft[norm_ucestanost3-4])+log(kmdfft[norm_ucestanost3+4]);}
		if(norm_ucestanost3>2)
		{energija_u_okolini_max3 = maksimum3+log(kmdfft[norm_ucestanost3-1])+log(kmdfft[norm_ucestanost3+1])+log(kmdfft[norm_ucestanost3-2])+log(kmdfft[norm_ucestanost3+2])+log(kmdfft[norm_ucestanost3-3])+log(kmdfft[norm_ucestanost3+3]);}
		if(norm_ucestanost3>1)
		{energija_u_okolini_max3 = maksimum3+log(kmdfft[norm_ucestanost3-1])+log(kmdfft[norm_ucestanost3+1])+log(kmdfft[norm_ucestanost3-2])+log(kmdfft[norm_ucestanost3+2]);}
		if(norm_ucestanost3==1)
		{energija_u_okolini_max3 = maksimum3+log(kmdfft[norm_ucestanost3-1])+log(kmdfft[norm_ucestanost3+1])+log(kmdfft[norm_ucestanost3+2]);}
		if(norm_ucestanost3==0)
		{energija_u_okolini_max3 = maksimum3+log(kmdfft[norm_ucestanost3+1])+log(kmdfft[norm_ucestanost3+2]);}
		vekt_obelezja[l][2]=energija_u_okolini_max3; log_energije[l][2]=energija_u_okolini_max3;
		}

		int norm_ucestanost4 = 1;
		if(broj_obelezja==4 || broj_obelezja==5 || broj_obelezja==6 || broj_obelezja==7 || broj_obelezja==8 || broj_obelezja==9 || broj_obelezja==10 || broj_obelezja==11 || broj_obelezja==12)
		{
		double maksimum4 = -10000; double energija_u_okolini_max4 = 0.0;
		for(int i=1;i<N/2;i++) //racunanje cetvrtog po redu maksimuma u spektru
		{
			if( ( ((i<norm_ucestanost1-5) || (i>norm_ucestanost1+5)) && ((i<norm_ucestanost2-5) || (i>norm_ucestanost2+5)) ) 
				&& ((i<norm_ucestanost3-5) || (i>norm_ucestanost3+5)) )
			{
				if(log(kmdfft[i])>maksimum4)
				{
					maksimum4 = log(kmdfft[i]); norm_ucestanost4=i;
				}
			}
		}
		vektor_ucestanosti_maksimuma4[l] = norm_ucestanost4;//ucestanost najvece komponente po cetvrtom kriterijumu u razmatranom prozoru
		vektor_ucestanosti_maksimuma[l][3]=norm_ucestanost4;
		if(norm_ucestanost4>1)
		{energija_u_okolini_max4 = maksimum4+log(kmdfft[norm_ucestanost4-1])+log(kmdfft[norm_ucestanost4+1])+log(kmdfft[norm_ucestanost4-2])+log(kmdfft[norm_ucestanost4+2]);}
		if(norm_ucestanost4==1)
		{energija_u_okolini_max4 = maksimum4+log(kmdfft[norm_ucestanost4-1])+log(kmdfft[norm_ucestanost4+1])+log(kmdfft[norm_ucestanost4+2]);}
		if(norm_ucestanost4==0)
		{energija_u_okolini_max4 = maksimum4+log(kmdfft[norm_ucestanost4+1])+log(kmdfft[norm_ucestanost4+2]);}
		vekt_obelezja[l][3]=energija_u_okolini_max4; log_energije[l][3]=energija_u_okolini_max4;
		}
		
		int norm_ucestanost5 = 1;
		if(broj_obelezja==5 || broj_obelezja==6 || broj_obelezja==7 || broj_obelezja==8 || broj_obelezja==9 || broj_obelezja==10 || broj_obelezja==11 || broj_obelezja==12)
		{
		double maksimum5 = -10000; double energija_u_okolini_max5 = 0.0;
		for(int i=1;i<N/2;i++) //racunanje petog po redu maksimuma u spektru
		{
			if( ( ((i<norm_ucestanost1-5) || (i>norm_ucestanost1+5)) && ((i<norm_ucestanost2-5) || (i>norm_ucestanost2+5)) ) 
				&& ((i<norm_ucestanost3-5) || (i>norm_ucestanost3+5)) && ((i<norm_ucestanost4-5) || (i>norm_ucestanost4+5)) )
			{
				if(log(kmdfft[i])>maksimum5)
				{
					maksimum5 = log(kmdfft[i]); norm_ucestanost5=i;
				}
			}
		}
		vektor_ucestanosti_maksimuma5[l] = norm_ucestanost5;//ucestanost najvece komponente po petom kriterijumu u razmatranom prozoru
		vektor_ucestanosti_maksimuma[l][4]=norm_ucestanost5;
		if(norm_ucestanost5>1)
		{energija_u_okolini_max5 = maksimum5+log(kmdfft[norm_ucestanost5-1])+log(kmdfft[norm_ucestanost5+1])+log(kmdfft[norm_ucestanost5-2])+log(kmdfft[norm_ucestanost5+2]);}
		if(norm_ucestanost5==1)
		{energija_u_okolini_max5 = maksimum5+log(kmdfft[norm_ucestanost5-1])+log(kmdfft[norm_ucestanost5+1])+log(kmdfft[norm_ucestanost5+2]);}
		if(norm_ucestanost5==0)
		{energija_u_okolini_max5 = maksimum5+log(kmdfft[norm_ucestanost5+1])+log(kmdfft[norm_ucestanost5+2]);}
		vekt_obelezja[l][4]=energija_u_okolini_max5; log_energije[l][4]=energija_u_okolini_max5;
		}

		int norm_ucestanost6 = 1;
		if(broj_obelezja==6 || broj_obelezja==7 || broj_obelezja==8 || broj_obelezja==9 || broj_obelezja==10 || broj_obelezja==11 || broj_obelezja==12)
		{
		double maksimum6 = -10000; double energija_u_okolini_max6 = 0.0;
		for(int i=1;i<N/2;i++) //racunanje sestog po redu maksimuma u spektru
		{
			if( ( ((i<norm_ucestanost1-5) || (i>norm_ucestanost1+5)) && ((i<norm_ucestanost2-5) || (i>norm_ucestanost2+5)) ) 
				&& ((i<norm_ucestanost3-5) || (i>norm_ucestanost3+5)) && ((i<norm_ucestanost4-5) || (i>norm_ucestanost4+5)) 
				&& ((i<norm_ucestanost5-5) || (i>norm_ucestanost5+5)) )
			{
				if(log(kmdfft[i])>maksimum6)
				{
					maksimum6 = log(kmdfft[i]); norm_ucestanost6=i;
				}
			}
		}
		vektor_ucestanosti_maksimuma6[l] = norm_ucestanost6;//ucestanost najvece komponente po sestom kriterijumu u razmatranom prozoru
		vektor_ucestanosti_maksimuma[l][5]=norm_ucestanost6;
		if(norm_ucestanost6>1)
		{energija_u_okolini_max6 = maksimum6+log(kmdfft[norm_ucestanost6-1])+log(kmdfft[norm_ucestanost6+1])+log(kmdfft[norm_ucestanost6-2])+log(kmdfft[norm_ucestanost6+2]);}
		if(norm_ucestanost6==1)
		{energija_u_okolini_max6 = maksimum6+log(kmdfft[norm_ucestanost6-1])+log(kmdfft[norm_ucestanost6+1])+log(kmdfft[norm_ucestanost6+2]);}
		if(norm_ucestanost6==0)
		{energija_u_okolini_max6 = maksimum6+log(kmdfft[norm_ucestanost6+1])+log(kmdfft[norm_ucestanost6+2]);}
		vekt_obelezja[l][5]=energija_u_okolini_max6; log_energije[l][5]=energija_u_okolini_max6;
		}

		int norm_ucestanost7 = 1;
		if(broj_obelezja==7 || broj_obelezja==8 || broj_obelezja==9 || broj_obelezja==10 || broj_obelezja==11 || broj_obelezja==12)
		{
		double maksimum7 = -10000; double energija_u_okolini_max7 = 0.0;
		for(int i=1;i<N/2;i++) //racunanje sedmog po redu maksimuma u spektru
		{
			if( ( ((i<norm_ucestanost1-5) || (i>norm_ucestanost1+5)) && ((i<norm_ucestanost2-5) || (i>norm_ucestanost2+5)) ) 
				&& ((i<norm_ucestanost3-5) || (i>norm_ucestanost3+5)) && ((i<norm_ucestanost4-5) || (i>norm_ucestanost4+5)) 
				&& ((i<norm_ucestanost5-5) || (i>norm_ucestanost5+5)) && ((i<norm_ucestanost6-5) || (i>norm_ucestanost6+5)) )
			{
				if(log(kmdfft[i])>maksimum7)
				{
					maksimum7 = log(kmdfft[i]); norm_ucestanost7=i;
				}
			}
		}
		vektor_ucestanosti_maksimuma7[l] = norm_ucestanost7;//ucestanost najvece komponente po sedmom kriterijumu u razmatranom prozoru
		vektor_ucestanosti_maksimuma[l][6]=norm_ucestanost7;
		if(norm_ucestanost7>1)
		{energija_u_okolini_max7 = maksimum7+log(kmdfft[norm_ucestanost7-1])+log(kmdfft[norm_ucestanost7+1])+log(kmdfft[norm_ucestanost7-2])+log(kmdfft[norm_ucestanost7+2]);}
		if(norm_ucestanost7==1)
		{energija_u_okolini_max7 = maksimum7+log(kmdfft[norm_ucestanost7-1])+log(kmdfft[norm_ucestanost7+1])+log(kmdfft[norm_ucestanost7+2]);}
		if(norm_ucestanost7==0)
		{energija_u_okolini_max7 = maksimum7+log(kmdfft[norm_ucestanost7+1])+log(kmdfft[norm_ucestanost7+2]);}
		vekt_obelezja[l][6]=energija_u_okolini_max7; log_energije[l][6]=energija_u_okolini_max7;
		}

		int norm_ucestanost8 = 1;
		if(broj_obelezja==8 || broj_obelezja==9 || broj_obelezja==10 || broj_obelezja==11 || broj_obelezja==12)
		{
		double maksimum8 = -10000; double energija_u_okolini_max8 = 0.0;
		for(int i=1;i<N/2;i++) //racunanje osmog po redu maksimuma u spektru
		{
			if( ( ((i<norm_ucestanost1-5) || (i>norm_ucestanost1+5)) && ((i<norm_ucestanost2-5) || (i>norm_ucestanost2+5)) ) 
				&& ((i<norm_ucestanost3-5) || (i>norm_ucestanost3+5)) && ((i<norm_ucestanost4-5) || (i>norm_ucestanost4+5)) 
				&& ((i<norm_ucestanost5-5) || (i>norm_ucestanost5+5)) && ((i<norm_ucestanost6-5) || (i>norm_ucestanost6+5)) 
				&& ((i<norm_ucestanost7-5) || (i>norm_ucestanost7+5)) )
			{
				if(log(kmdfft[i])>maksimum8)
				{
					maksimum8 = log(kmdfft[i]); norm_ucestanost8=i;
				}
			}
		}
		vektor_ucestanosti_maksimuma8[l] = norm_ucestanost8;//ucestanost najvece komponente po osmom kriterijumu u razmatranom prozoru
		vektor_ucestanosti_maksimuma[l][7]=norm_ucestanost8;
		if(norm_ucestanost8>1)
		{energija_u_okolini_max8 = maksimum8+log(kmdfft[norm_ucestanost8-1])+log(kmdfft[norm_ucestanost8+1])+log(kmdfft[norm_ucestanost8-2])+log(kmdfft[norm_ucestanost8+2]);}
		if(norm_ucestanost8==1)
		{energija_u_okolini_max8 = maksimum8+log(kmdfft[norm_ucestanost8-1])+log(kmdfft[norm_ucestanost8+1])+log(kmdfft[norm_ucestanost8+2]);}
		if(norm_ucestanost8==0)
		{energija_u_okolini_max8 = maksimum8+log(kmdfft[norm_ucestanost8+1])+log(kmdfft[norm_ucestanost8+2]);}
		vekt_obelezja[l][7]=energija_u_okolini_max8; log_energije[l][7]=energija_u_okolini_max8;
		}

		int norm_ucestanost9 = 1;
		if(broj_obelezja==9 || broj_obelezja==10 || broj_obelezja==11 || broj_obelezja==12)
		{
		double maksimum9 = -10000; double energija_u_okolini_max9 = 0.0;
		for(int i=1;i<N/2;i++) //racunanje devetog po redu maksimuma u spektru
		{
			if( ( ((i<norm_ucestanost1-5) || (i>norm_ucestanost1+5)) && ((i<norm_ucestanost2-5) || (i>norm_ucestanost2+5)) ) 
				&& ((i<norm_ucestanost3-5) || (i>norm_ucestanost3+5)) && ((i<norm_ucestanost4-5) || (i>norm_ucestanost4+5)) 
				&& ((i<norm_ucestanost5-5) || (i>norm_ucestanost5+5)) && ((i<norm_ucestanost6-5) || (i>norm_ucestanost6+5)) 
				&& ((i<norm_ucestanost7-5) || (i>norm_ucestanost7+5)) && ((i<norm_ucestanost8-5) || (i>norm_ucestanost8+5)) )
			{
				if(log(kmdfft[i])>maksimum9)
				{
					maksimum9 = log(kmdfft[i]); norm_ucestanost9=i;
				}
			}
		}
		vektor_ucestanosti_maksimuma9[l] = norm_ucestanost9;//ucestanost najvece komponente po devetom kriterijumu u razmatranom prozoru
		vektor_ucestanosti_maksimuma[l][8]=norm_ucestanost9;
		if(norm_ucestanost9>1)
		{energija_u_okolini_max9 = maksimum9+log(kmdfft[norm_ucestanost9-1])+log(kmdfft[norm_ucestanost9+1])+log(kmdfft[norm_ucestanost9-2])+log(kmdfft[norm_ucestanost9+2]);}
		if(norm_ucestanost9==1)
		{energija_u_okolini_max9 = maksimum9+log(kmdfft[norm_ucestanost9-1])+log(kmdfft[norm_ucestanost9+1])+log(kmdfft[norm_ucestanost9+2]);}
		if(norm_ucestanost9==0)
		{energija_u_okolini_max9 = maksimum9+log(kmdfft[norm_ucestanost9+1])+log(kmdfft[norm_ucestanost9+2]);}
		vekt_obelezja[l][8]=energija_u_okolini_max9; log_energije[l][8]=energija_u_okolini_max9;
		}

		int norm_ucestanost10 = 1;
		if(broj_obelezja==10 || broj_obelezja==11 || broj_obelezja==12)
		{
		double maksimum10 = -10000; double energija_u_okolini_max10 = 0.0;
		for(int i=1;i<N/2;i++) //racunanje desetog po redu maksimuma u spektru
		{
			if( ( ((i<norm_ucestanost1-5) || (i>norm_ucestanost1+5)) && ((i<norm_ucestanost2-5) || (i>norm_ucestanost2+5)) ) 
				&& ((i<norm_ucestanost3-5) || (i>norm_ucestanost3+5)) && ((i<norm_ucestanost4-5) || (i>norm_ucestanost4+5)) 
				&& ((i<norm_ucestanost5-5) || (i>norm_ucestanost5+5)) && ((i<norm_ucestanost6-5) || (i>norm_ucestanost6+5)) 
				&& ((i<norm_ucestanost7-5) || (i>norm_ucestanost7+5)) && ((i<norm_ucestanost8-5) || (i>norm_ucestanost8+5)) 
				&& ((i<norm_ucestanost9-5) || (i>norm_ucestanost9+5)) )
			{
				if(log(kmdfft[i])>maksimum10)
				{
					maksimum10 = log(kmdfft[i]); norm_ucestanost10=i;
				}
			}
		}
		vektor_ucestanosti_maksimuma10[l] = norm_ucestanost10;//ucestanost najvece komponente po desetom kriterijumu u razmatranom prozoru
		vektor_ucestanosti_maksimuma[l][9]=norm_ucestanost10;
		if(norm_ucestanost10>1)
		{energija_u_okolini_max10 = maksimum10+log(kmdfft[norm_ucestanost10-1])+log(kmdfft[norm_ucestanost10+1])+log(kmdfft[norm_ucestanost10-2])+log(kmdfft[norm_ucestanost10+2]);}
		if(norm_ucestanost10==1)
		{energija_u_okolini_max10 = maksimum10+log(kmdfft[norm_ucestanost10-1])+log(kmdfft[norm_ucestanost10+1])+log(kmdfft[norm_ucestanost10+2]);}
		if(norm_ucestanost10==0)
		{energija_u_okolini_max10 = maksimum10+log(kmdfft[norm_ucestanost10+1])+log(kmdfft[norm_ucestanost10+2]);}
		vekt_obelezja[l][9]=energija_u_okolini_max10; log_energije[l][9]=energija_u_okolini_max10;
		}

		int norm_ucestanost11 = 1;
		if(broj_obelezja==11 || broj_obelezja==12)
		{
		double maksimum11 = -10000; double energija_u_okolini_max11 = 0.0;
		for(int i=1;i<N/2;i++) //racunanje 11. po redu maksimuma u spektru
		{
			if( ( ((i<norm_ucestanost1-5) || (i>norm_ucestanost1+5)) && ((i<norm_ucestanost2-5) || (i>norm_ucestanost2+5)) ) 
				&& ((i<norm_ucestanost3-5) || (i>norm_ucestanost3+5)) && ((i<norm_ucestanost4-5) || (i>norm_ucestanost4+5)) 
				&& ((i<norm_ucestanost5-5) || (i>norm_ucestanost5+5)) && ((i<norm_ucestanost6-5) || (i>norm_ucestanost6+5)) 
				&& ((i<norm_ucestanost7-5) || (i>norm_ucestanost7+5)) && ((i<norm_ucestanost8-5) || (i>norm_ucestanost8+5)) 
				&& ((i<norm_ucestanost9-5) || (i>norm_ucestanost9+5)) && ((i<norm_ucestanost10-5)||(i>norm_ucestanost10+5)) )
			{
				if(log(kmdfft[i])>maksimum11)
				{
					maksimum11 = log(kmdfft[i]); norm_ucestanost11=i;
				}
			}
		}
		vektor_ucestanosti_maksimuma11[l] = norm_ucestanost11;//ucestanost najvece komponente po 11. kriterijumu u razmatranom prozoru
		vektor_ucestanosti_maksimuma[l][10]=norm_ucestanost11;
		if(norm_ucestanost11>1)
		{energija_u_okolini_max11 = maksimum11+log(kmdfft[norm_ucestanost11-1])+log(kmdfft[norm_ucestanost11+1])+log(kmdfft[norm_ucestanost11-2])+log(kmdfft[norm_ucestanost11+2]);}
		if(norm_ucestanost11==1)
		{energija_u_okolini_max11 = maksimum11+log(kmdfft[norm_ucestanost11-1])+log(kmdfft[norm_ucestanost11+1])+log(kmdfft[norm_ucestanost11+2]);}
		if(norm_ucestanost11==0)
		{energija_u_okolini_max11 = maksimum11+log(kmdfft[norm_ucestanost11+1])+log(kmdfft[norm_ucestanost11+2]);}
		vekt_obelezja[l][10]=energija_u_okolini_max11; log_energije[l][10]=energija_u_okolini_max11;
		}

		int norm_ucestanost12 = 1;
		if(broj_obelezja==12)
		{
		double maksimum12 = -10000; double energija_u_okolini_max12 = 0.0;
		for(int i=1;i<N/2;i++) //racunanje 12. po redu maksimuma u spektru
		{
			if( ( ((i<norm_ucestanost1-5) || (i>norm_ucestanost1+5)) && ((i<norm_ucestanost2-5) || (i>norm_ucestanost2+5)) ) 
				&& ((i<norm_ucestanost3-5) || (i>norm_ucestanost3+5)) && ((i<norm_ucestanost4-5) || (i>norm_ucestanost4+5)) 
				&& ((i<norm_ucestanost5-5) || (i>norm_ucestanost5+5)) && ((i<norm_ucestanost6-5) || (i>norm_ucestanost6+5)) 
				&& ((i<norm_ucestanost7-5) || (i>norm_ucestanost7+5)) && ((i<norm_ucestanost8-5) || (i>norm_ucestanost8+5)) 
				&& ((i<norm_ucestanost9-5) || (i>norm_ucestanost9+5)) && ((i<norm_ucestanost10-5)||(i>norm_ucestanost10+5)) 
				&& ((i<norm_ucestanost11-5) || (i>norm_ucestanost11+5)) )
			{
				if(log(kmdfft[i])>maksimum12)
				{
					maksimum12 = log(kmdfft[i]); norm_ucestanost12=i;
				}
			}
		}
		vektor_ucestanosti_maksimuma12[l] = norm_ucestanost12;//ucestanost najvece komponente po 12. kriterijumu u razmatranom prozoru
		vektor_ucestanosti_maksimuma[l][11]=norm_ucestanost12;
		if(norm_ucestanost12>1)
		{energija_u_okolini_max12 = maksimum12+log(kmdfft[norm_ucestanost12-1])+log(kmdfft[norm_ucestanost12+1])+log(kmdfft[norm_ucestanost12-2])+log(kmdfft[norm_ucestanost12+2]);}
		if(norm_ucestanost12==1)
		{energija_u_okolini_max12 = maksimum12+log(kmdfft[norm_ucestanost12-1])+log(kmdfft[norm_ucestanost12+1])+log(kmdfft[norm_ucestanost12+2]);}
		if(norm_ucestanost12==0)
		{energija_u_okolini_max12 = maksimum12+log(kmdfft[norm_ucestanost12+1])+log(kmdfft[norm_ucestanost12+2]);}
		vekt_obelezja[l][11]=energija_u_okolini_max12; log_energije[l][11]=energija_u_okolini_max12;
		}
	}//l

	//ponderisanje dodatnih obelezja shodno njihovim vrednostima u odnosu na najveci maksimum po prozorima
	for(int i=0;i<brojprozora;i++)
	{
		if(broj_obelezja==1 || broj_obelezja==2 || broj_obelezja==3 || broj_obelezja==4 || broj_obelezja==5)
		{vekt_obelezja[i][0] = (vekt_obelezja[i][0]*vekt_obelezja[i][0]/najveci_maksimum_po_prozorima)/((double)vektor_ucestanosti_maksimuma1[i]);}
		if(broj_obelezja==2 || broj_obelezja==3 || broj_obelezja==4 || broj_obelezja==5)
		{vekt_obelezja[i][1] = (vekt_obelezja[i][1]*vekt_obelezja[i][1]/najveci_maksimum_po_prozorima)/((double)vektor_ucestanosti_maksimuma2[i]);}	
		if(broj_obelezja==3 || broj_obelezja==4 || broj_obelezja==5)
		{vekt_obelezja[i][2] = (vekt_obelezja[i][2]*vekt_obelezja[i][2]/najveci_maksimum_po_prozorima)/((double)vektor_ucestanosti_maksimuma3[i]);}
		if(broj_obelezja==4 || broj_obelezja==5)
		{vekt_obelezja[i][3] = (vekt_obelezja[i][3]*vekt_obelezja[i][3]/najveci_maksimum_po_prozorima)/((double)vektor_ucestanosti_maksimuma4[i]);}
		if(broj_obelezja==5)
		{vekt_obelezja[i][4] = (vekt_obelezja[i][4]*vekt_obelezja[i][4]/najveci_maksimum_po_prozorima)/((double)vektor_ucestanosti_maksimuma5[i]);}
		if(broj_obelezja==6 || broj_obelezja==7 || broj_obelezja==8 || broj_obelezja==9 || broj_obelezja==10 || broj_obelezja==11 || broj_obelezja==12)
		{vekt_obelezja[i][5] = (vekt_obelezja[i][5]*vekt_obelezja[i][5]/najveci_maksimum_po_prozorima)/((double)vektor_ucestanosti_maksimuma6[i]);}
		//if(broj_obelezja==7)
		//{vekt_obelezja[i][6] = (vekt_obelezja[i][6]*vekt_obelezja[i][6]/najveci_maksimum_po_prozorima)/((double)vektor_ucestanosti_maksimuma7[i]);}
		//if(broj_obelezja==8)
		//{vekt_obelezja[i][7] = (vekt_obelezja[i][7]*vekt_obelezja[i][7]/najveci_maksimum_po_prozorima)/((double)vektor_ucestanosti_maksimuma8[i]);}
		for(int j=0;j<broj_obelezja;j++)
		{vekt_obelezja[i][j]=0.0; 
		for(int k=0;k<10;k++)
		//{vekt_obelezja[i][j]+=log_energije[i][k]*cos(((double)j)* (((double)(k+1))-0.5));}
		//{vekt_obelezja[i][j]+=log_energije[i][k]*cos(((double)j)* (((double)(vektor_ucestanosti_maksimuma[i][k]))-0.5));}
		{vekt_obelezja[i][j]+=log_energije[i][k]*cos(((double)j)* (((double)(vektor_ucestanosti_maksimuma[i][k]+1))-0.5)*PI/(double)vektor_ucestanosti_maksimuma[i][9]);} 
		}
	}
	return vekt_obelezja;
}

vector<vector<double>> dmfcc(vector<vector<double>> vektoriobelezja, int brojMFCCs)
{
	vector<vector<double>> dvektoriobelezja(brojprozora, vector<double> (brojMFCCs));
	for(int i=0; i<brojprozora; i++)
	{
		for(int j=0; j<brojMFCCs; j++)
		{
			dvektoriobelezja[0][j]=0;
			dvektoriobelezja[1][j]=0;
			if(i>=2 && i<=(brojprozora-3))
			{
				dvektoriobelezja[i][j]=(vektoriobelezja[i][j]+vektoriobelezja[i+1][j]+2*vektoriobelezja[i+2][j]-vektoriobelezja[i-1][j]-2*vektoriobelezja[i-2][j])/6;
			}
			dvektoriobelezja[brojprozora-2][j]=0;
			dvektoriobelezja[brojprozora-1][j]=0;
		}
	}
	return dvektoriobelezja;
}

vector<vector<double>> ddmfcc(vector<vector<double>> dvektoriobelezja,int brojMFCCs)
{
	vector<vector<double>> ddvektoriobelezja(brojprozora, vector<double> (brojMFCCs));
	for(int i=0; i<brojprozora; i++)
	{
		for(int j=0; j<brojMFCCs; j++)
		{
			ddvektoriobelezja[0][j]=0;
			ddvektoriobelezja[1][j]=0;
			if(i>=2 && i<=(brojprozora-3))
			{
				ddvektoriobelezja[i][j]=(dvektoriobelezja[i][j]+dvektoriobelezja[i+1][j]+2*dvektoriobelezja[i+2][j]-dvektoriobelezja[i-1][j]-2*dvektoriobelezja[i-2][j])/6;
			}
			ddvektoriobelezja[brojprozora-2][j]=0;
			ddvektoriobelezja[brojprozora-1][j]=0;
		}
	}
	return ddvektoriobelezja;
}

vector<vector<double>> dddmfcc(vector<vector<double>> ddvektoriobelezja,int brojMFCCs)
{
	vector<vector<double>> dddvektoriobelezja(brojprozora, vector<double> (brojMFCCs));
	for(int i=0; i<brojprozora; i++)
	{
		for(int j=0; j<brojMFCCs; j++)
		{
			dddvektoriobelezja[0][j]=0;
			dddvektoriobelezja[1][j]=0;
			if(i>=2 && i<=(brojprozora-3))
			{
				dddvektoriobelezja[i][j]=(ddvektoriobelezja[i][j]+ddvektoriobelezja[i+1][j]+2*ddvektoriobelezja[i+2][j]-ddvektoriobelezja[i-1][j]-2*ddvektoriobelezja[i-2][j])/6;
			}
			dddvektoriobelezja[brojprozora-2][j]=0;
			dddvektoriobelezja[brojprozora-1][j]=0;
		}
	}
	return dddvektoriobelezja;
}

vector<vector<double>> ddddmfcc(vector<vector<double>> dddvektoriobelezja,int brojMFCCs)
{
	vector<vector<double>> ddddvektoriobelezja(brojprozora, vector<double> (brojMFCCs));
	for(int i=0; i<brojprozora; i++)
	{
		for(int j=0; j<brojMFCCs; j++)
		{
			ddddvektoriobelezja[0][j]=0;
			ddddvektoriobelezja[1][j]=0;
			if(i>=2 && i<=(brojprozora-3))
			{
				ddddvektoriobelezja[i][j]=(dddvektoriobelezja[i][j]+dddvektoriobelezja[i+1][j]+2*dddvektoriobelezja[i+2][j]-dddvektoriobelezja[i-1][j]-2*dddvektoriobelezja[i-2][j])/6;
			}
			ddddvektoriobelezja[brojprozora-2][j]=0;
			ddddvektoriobelezja[brojprozora-1][j]=0;
		}
	}
	return ddddvektoriobelezja;
}


vector<double> veksrvr(int brojMFCCs,vector<vector<double>> dvektoriobelezja)
{
	vector<double> veksrvr(brojMFCCs);
	for (int j=0; j<brojMFCCs; j++)
	{
		for (int i=0; i<brojprozora; i++)
		{
			veksrvr[j] += (dvektoriobelezja[i] [j])/brojprozora;
		}
	}
	return veksrvr;
}

vector<vector<double>> kovmat(int brojMFCCs,vector<vector<double> > dvektoriobelezja,vector<double> veksrvr)
{
	vector<vector<double>> kovmat(brojMFCCs, vector<double> (brojMFCCs)); //kovarijansna matrica
	for (int i=0; i<brojMFCCs; i++)
	{
		for (int j=0; j<brojMFCCs; j++)
		{
			for (int k=0; k<brojprozora; k++)
			{
				kovmat[i][j] += (((dvektoriobelezja[k][i]-veksrvr[i]) * (dvektoriobelezja[k][j]-veksrvr[j]))/(brojprozora-1));
			}
		}
	}
	return kovmat;
}

vector<double> ucitavanjefajlova(int brojfajlova, vector<const char*> fajlovi)
{
	vector<double> odbirci;
	vector<double> priv;
	odbirci.swap(ucitavanjewavfajla(fajlovi[0]));
	for(int i=1; i<brojfajlova; i++)
	{
		priv.swap(ucitavanjewavfajla(fajlovi[i]));
		odbirci.insert(odbirci.end(),priv.begin(),priv.end());
	}
	return odbirci;
}


double apsodstupanjesrvr(int brojMFCCs, vector<double> veksrvr, vector<double> veksrvrt)
{
	double apsodstsrvr=0;
	for (int i=0; i<brojMFCCs; i++)
	{
		apsodstsrvr += (abs(veksrvr[i] - veksrvrt[i]))/12;
	}
	return apsodstsrvr;
}
double apsodstupanjekovmat(int brojMFCCs, vector<vector<double>> kovmat, vector<vector<double>> kovmatt)
{
	double apsodstupanjekovmatrice=0;
	for(int i=0; i<brojMFCCs; i++)  //elementi kovarijansne matrice su simetricni oko glavne dijagonale
	{            //stoga dovoljno je posmatrati razlikovanje elemenata na glavnoj dijagonali i elemenata ispod ili iznad nje  
		for(int j=0; j<=i; j++)//j<i+1 ili moze j<=i umesto j<brojMFCCs
		{
			if(j==i){
				apsodstupanjekovmatrice += (abs(kovmat[i][j] - kovmatt[i][j]))/(brojMFCCs*brojMFCCs);}
			else {apsodstupanjekovmatrice += ((abs(kovmat[i][j] - kovmatt[i][j]))*2.0)/(brojMFCCs*brojMFCCs);}
		}//za j<i+1 ili j<=i posmatrali bi se elementi ispod glavne dijagonale i na njoj
	}
	return apsodstupanjekovmatrice;
}

int odlucivanje(int brojgovornika, vector<string> vekidgovornika, vector<vector<double>> vekapsodstsrvr, vector<vector<double>> vekapsodstkovmat)
{
	int prepoznati;  //redni broj prepoznatog govornika
	int brojtacnoprepoznatih = 0;
	float tacnost;
	int prepoznati1; //redni broj drugog najslicnijeg govornika posmatranom govorniku
	for(int i=0;i<brojgovornika;i++)  //test govornik (prepoznavani glas)
	{
		float min=vekapsodstkovmat[i][0];
		float min1=vekapsodstkovmat[i][0];
		prepoznati=0;
		prepoznati1=0;
		for (int j=1; j<brojgovornika; j++)   // redni broj pretpostavljenog govornika, govornika ciji su modeli dobijeni prilikom obuke
		{
			if (vekapsodstkovmat[i][j]<min)
			{
				min1=min; //radi cuvanja prvog veceg razlikovanja u odnosu na ustanovljeno minimalno
				min=vekapsodstkovmat[i][j];
				prepoznati1=prepoznati; //drugi najslicniji govornik posmatranom govorniku
				prepoznati=j;
			}
		}
		cout << "\nGovornik " << vekidgovornika[i] << " prepoznat je kao govornik " << vekidgovornika[prepoznati] << " pri cemu mera razlikovanja iznosi " << min << " a prvo vece mera razlikovanja iznosi " << min1 << " u odnosu na govornika " << vekidgovornika[prepoznati1] <<  endl;
		if (prepoznati==i)
		{
			brojtacnoprepoznatih += 1;
		}
	}
	tacnost = (float)brojtacnoprepoznatih/(float)brojgovornika; //-(float)broj_nepouzdanih;
	cout << "Tacnost prepoznavanja iznosi: " << brojtacnoprepoznatih << "/" << brojgovornika << " = " << tacnost*100.0 << "%" << endl;
	return 0;
}

int prepoznavanje_nad_snimkom(char* snimak, int brojgovornika, vector<string> vekidgovornika, vector<double> vekapsodstkovmat_test)
{
	int prepoznati=0,prepoznati1=0;   //pocetne vrednosti rednih brojeva: prepoznatog govornika i drugog najslicnijeg govornika posmatranom govorniku
	double min=vekapsodstkovmat_test[0]; double min1=vekapsodstkovmat_test[0];
	for(int i=0;i<brojgovornika;i++)  //radi poredjenja sa modelima svih govornika u govornoj bazi
	{
		if (vekapsodstkovmat_test[i]<min)
		{
			min1=min; //radi cuvanja prvog veceg razlikovanja u odnosu na ustanovljeno minimalno
			min=vekapsodstkovmat_test[i];
			prepoznati1=prepoznati; //drugi najslicniji govornik posmatranom govorniku
			prepoznati=i;
		}
	}
	cout << "\nSnimak " << snimak << " prepoznat je kao govornik " << vekidgovornika[prepoznati] << " pri cemu mera razlikovanja iznosi " << min << " a prva veca mera razlikovanja iznosi " << min1 << " u odnosu na govornika " << vekidgovornika[prepoznati1] <<  endl;
	//if(vekidgovornika[prepoznati]=="StevanJokic" | vekidgovornika[prepoznati]=="IvanJokic")
	//{system("\"C:\\Program Files\\Winamp\\winamp.exe\"/NEW Vo_Kosi_Da_Ti_Spijam.mp3");}
		return 0;
}

string prepoznavanje_nad_snimkom_s(string snimak, int brojgovornika, vector<string> vekidgovornika, vector<double> vekapsodstkovmat_test)
{
	string prepoznati_s = "";
	int prepoznati=0,prepoznati1=0;   //pocetne vrednosti rednih brojeva: prepoznatog govornika i drugog najslicnijeg govornika posmatranom govorniku
	double min=vekapsodstkovmat_test[0]; double min1=vekapsodstkovmat_test[0];
	for(int i=0;i<brojgovornika;i++)  //radi poredjenja sa modelima svih govornika u govornoj bazi
	{
		if (vekapsodstkovmat_test[i]<min)
		{
			min1=min; //radi cuvanja prvog veceg razlikovanja u odnosu na ustanovljeno minimalno
			min=vekapsodstkovmat_test[i];
			prepoznati1=prepoznati; //drugi najslicniji govornik posmatranom govorniku
			prepoznati=i;
		}
	}
	cout << "\nSnimak " << snimak << " prepoznat je kao govornik " << vekidgovornika[prepoznati] << " pri cemu mera razlikovanja iznosi " << min << " a prva veca mera razlikovanja iznosi " << min1 << " u odnosu na govornika " << vekidgovornika[prepoznati1] <<  endl;
	//if(vekidgovornika[prepoznati]=="StevanJokic" | vekidgovornika[prepoznati]=="IvanJokic")
	//{system("\"C:\\Program Files\\Winamp\\winamp.exe\"/NEW Vo_Kosi_Da_Ti_Spijam.mp3");}
	prepoznati_s = vekidgovornika[prepoznati];	
	return prepoznati_s;
}

int knn_prepoznavanje_nad_snimkom(char* snimak, int brojgovornika, vector<string> vekidgovornika, 
								  vector<double> vekapsodstkovmat_test)
{
	int prepoznati1=0; double min1=vekapsodstkovmat_test[0];
	for(int i=0;i<brojgovornika;i++)  //radi poredjenja sa modelima svih govornika u govornoj bazi
	{
		if (vekapsodstkovmat_test[i]<min1)
		{
			min1=vekapsodstkovmat_test[i];
			prepoznati1=i;
		}
	}
	vekapsodstkovmat_test[prepoznati1]=1000.0; //na mesto prethodno odredjenog najmanjeg rastojanja
	          //postavlja se izuzetno velika vrednost za koju se zna da ce biti najveca i zato nece 
		      // uticati na odluku o sledecem najmanjem elementu
	
	int prepoznati2=0; double min2 = vekapsodstkovmat_test[0];
	for(int i=0;i<brojgovornika;i++)  //radi poredjenja sa modelima svih govornika u govornoj bazi
	{
		if (vekapsodstkovmat_test[i]<min2)
		{
			min2=vekapsodstkovmat_test[i];
			prepoznati2=i;
		}
	}
	vekapsodstkovmat_test[prepoznati2]=1000.0;
	
	int prepoznati3=0; double min3 = vekapsodstkovmat_test[0];
	for(int i=0;i<brojgovornika;i++)  //radi poredjenja sa modelima svih govornika u govornoj bazi
	{
		if (vekapsodstkovmat_test[i]<min3)
		{
			min3=vekapsodstkovmat_test[i];
			prepoznati3=i;
		}
	}
	vekapsodstkovmat_test[prepoznati3]=1000.0;
	
	int prepoznati4=0; double min4 = vekapsodstkovmat_test[0];
	for(int i=0;i<brojgovornika;i++)  //radi poredjenja sa modelima svih govornika u govornoj bazi
	{
		if (vekapsodstkovmat_test[i]<min4)
		{
			min4=vekapsodstkovmat_test[i];
			prepoznati4=i;
		}
	}
	vekapsodstkovmat_test[prepoznati4]=1000.0;
	
	int prepoznati5=0; double min5 = vekapsodstkovmat_test[0];
	for(int i=0;i<brojgovornika;i++)  //radi poredjenja sa modelima svih govornika u govornoj bazi
	{
		if (vekapsodstkovmat_test[i]<min5)
		{
			min5=vekapsodstkovmat_test[i];
			prepoznati5=i;
		}
	}
	vekapsodstkovmat_test[prepoznati5]=1000.0;
	
	int prepoznati6=0; double min6 = vekapsodstkovmat_test[0];
	for(int i=0;i<brojgovornika;i++)  //radi poredjenja sa modelima svih govornika u govornoj bazi
	{
		if (vekapsodstkovmat_test[i]<min6)
		{
			min6=vekapsodstkovmat_test[i];
			prepoznati6=i;
		}
	}
	vekapsodstkovmat_test[prepoznati6]=1000.0;
	
	int prepoznati7=0; double min7 = vekapsodstkovmat_test[0];
	for(int i=0;i<brojgovornika;i++)  //radi poredjenja sa modelima svih govornika u govornoj bazi
	{
		if (vekapsodstkovmat_test[i]<min7)
		{
			min7=vekapsodstkovmat_test[i];
			prepoznati7=i;
		}
	}
	vekapsodstkovmat_test[prepoznati7]=1000.0;
	
	int prepoznati8=0; double min8 = vekapsodstkovmat_test[0];
	for(int i=0;i<brojgovornika;i++)  //radi poredjenja sa modelima svih govornika u govornoj bazi
	{
		if (vekapsodstkovmat_test[i]<min8)
		{
			min8=vekapsodstkovmat_test[i];
			prepoznati8=i;
		}
	}
	vekapsodstkovmat_test[prepoznati8]=1000.0;
	
	int prepoznati9=0; double min9 = vekapsodstkovmat_test[0];
	for(int i=0;i<brojgovornika;i++)  //radi poredjenja sa modelima svih govornika u govornoj bazi
	{
		if (vekapsodstkovmat_test[i]<min9)
		{
			min9=vekapsodstkovmat_test[i];
			prepoznati9=i;
		}
	}
	vekapsodstkovmat_test[prepoznati9]=1000.0;
	
	int prepoznati10=0; double min10 = vekapsodstkovmat_test[0];
	for(int i=0;i<brojgovornika;i++)  //radi poredjenja sa modelima svih govornika u govornoj bazi
	{
		if (vekapsodstkovmat_test[i]<min10)
		{
			min10=vekapsodstkovmat_test[i];
			prepoznati10=i;
		}
	}
	/*vekapsodstkovmat_test[prepoznati10]=1000.0;
	
	int prepoznati11=0; double min11 = vekapsodstkovmat_test[0];
	for(int i=0;i<brojgovornika;i++)  //radi poredjenja sa modelima svih govornika u govornoj bazi
	{
		if (vekapsodstkovmat_test[i]<min11)
		{
			min11=vekapsodstkovmat_test[i];
			prepoznati11=i;
		}
	}
	vekapsodstkovmat_test[prepoznati11]=1000.0;
	
	int prepoznati12=0; double min12 = vekapsodstkovmat_test[0];
	for(int i=0;i<brojgovornika;i++)  //radi poredjenja sa modelima svih govornika u govornoj bazi
	{
		if (vekapsodstkovmat_test[i]<min12)
		{
			min12=vekapsodstkovmat_test[i];
			prepoznati12=i;
		}
	}
	vekapsodstkovmat_test[prepoznati12]=1000.0;
	
	int prepoznati13=0; double min13 = vekapsodstkovmat_test[0];
	for(int i=0;i<brojgovornika;i++)  //radi poredjenja sa modelima svih govornika u govornoj bazi
	{
		if (vekapsodstkovmat_test[i]<min13)
		{
			min13=vekapsodstkovmat_test[i];
			prepoznati13=i;
		}
	}
	vekapsodstkovmat_test[prepoznati13]=1000.0;
	
	int prepoznati14=0; double min14 = vekapsodstkovmat_test[0];
	for(int i=0;i<brojgovornika;i++)  //radi poredjenja sa modelima svih govornika u govornoj bazi
	{
		if (vekapsodstkovmat_test[i]<min14)
		{
			min14=vekapsodstkovmat_test[i];
			prepoznati14=i;
		}
	}
	vekapsodstkovmat_test[prepoznati14]=1000.0;
	
	int prepoznati15=0; double min15 = vekapsodstkovmat_test[0];
	for(int i=0;i<brojgovornika;i++)  //radi poredjenja sa modelima svih govornika u govornoj bazi
	{
		if (vekapsodstkovmat_test[i]<min15)
		{
			min15=vekapsodstkovmat_test[i];
			prepoznati15=i;
		}
	}
	vekapsodstkovmat_test[prepoznati15]=1000.0;
	
	int prepoznati16=0; double min16 = vekapsodstkovmat_test[0];
	for(int i=0;i<brojgovornika;i++)  //radi poredjenja sa modelima svih govornika u govornoj bazi
	{
		if (vekapsodstkovmat_test[i]<min16)
		{
			min16=vekapsodstkovmat_test[i];
			prepoznati16=i;
		}
	}
	vekapsodstkovmat_test[prepoznati16]=1000.0;
	
	int prepoznati17=0; double min17 = vekapsodstkovmat_test[0];
	for(int i=0;i<brojgovornika;i++)  //radi poredjenja sa modelima svih govornika u govornoj bazi
	{
		if (vekapsodstkovmat_test[i]<min17)
		{
			min17=vekapsodstkovmat_test[i];
			prepoznati17=i;
		}
	}
	vekapsodstkovmat_test[prepoznati17]=1000.0;
	
	int prepoznati18=0; double min18 = vekapsodstkovmat_test[0];
	for(int i=0;i<brojgovornika;i++)  //radi poredjenja sa modelima svih govornika u govornoj bazi
	{
		if (vekapsodstkovmat_test[i]<min18)
		{
			min18=vekapsodstkovmat_test[i];
			prepoznati18=i;
		}
	}
	vekapsodstkovmat_test[prepoznati18]=1000.0;
	
	int prepoznati19=0; double min19 = vekapsodstkovmat_test[0];
	for(int i=0;i<brojgovornika;i++)  //radi poredjenja sa modelima svih govornika u govornoj bazi
	{
		if (vekapsodstkovmat_test[i]<min19)
		{
			min19=vekapsodstkovmat_test[i];
			prepoznati19=i;
		}
	}
	vekapsodstkovmat_test[prepoznati19]=1000.0;
	
	int prepoznati20=0; double min20 = vekapsodstkovmat_test[0];
	for(int i=0;i<brojgovornika;i++)  //radi poredjenja sa modelima svih govornika u govornoj bazi
	{
		if (vekapsodstkovmat_test[i]<min20)
		{
			min20=vekapsodstkovmat_test[i];
			prepoznati20=i;
		}
	}*/
	
	cout << "\nZa snimak " << snimak << endl 
		<< " 1. najslicniji: " << vekidgovornika[prepoznati1] << " pri cemu mera razlikovanja iznosi " << min1 << endl
		<< " 2. najslicniji: " << vekidgovornika[prepoznati2] << " pri cemu mera razlikovanja iznosi " << min2 << endl
		<< " 3. najslicniji: " << vekidgovornika[prepoznati3] << " pri cemu mera razlikovanja iznosi " << min3 << endl
		<< " 4. najslicniji: " << vekidgovornika[prepoznati4] << " pri cemu mera razlikovanja iznosi " << min4 << endl
		<< " 5. najslicniji: " << vekidgovornika[prepoznati5] << " pri cemu mera razlikovanja iznosi " << min5 << endl
		<< " 6. najslicniji: " << vekidgovornika[prepoznati6] << " pri cemu mera razlikovanja iznosi " << min6 << endl
		<< " 7. najslicniji: " << vekidgovornika[prepoznati7] << " pri cemu mera razlikovanja iznosi " << min7 << endl
		<< " 8. najslicniji: " << vekidgovornika[prepoznati8] << " pri cemu mera razlikovanja iznosi " << min8 << endl
		<< " 9. najslicniji: " << vekidgovornika[prepoznati9] << " pri cemu mera razlikovanja iznosi " << min9 << endl
		<< " 10. najslicniji: " << vekidgovornika[prepoznati10] << " pri cemu mera razlikovanja iznosi " << min10 << endl;
		/*<< " 11. najslicniji: " << vekidgovornika[prepoznati11] << " pri cemu mera razlikovanja iznosi " << min11 << endl
		<< " 12. najslicniji: " << vekidgovornika[prepoznati12] << " pri cemu mera razlikovanja iznosi " << min12 << endl
		<< " 13. najslicniji: " << vekidgovornika[prepoznati13] << " pri cemu mera razlikovanja iznosi " << min13 << endl
		<< " 14. najslicniji: " << vekidgovornika[prepoznati14] << " pri cemu mera razlikovanja iznosi " << min14 << endl
		<< " 15. najslicniji: " << vekidgovornika[prepoznati15] << " pri cemu mera razlikovanja iznosi " << min15 << endl
		<< " 16. najslicniji: " << vekidgovornika[prepoznati16] << " pri cemu mera razlikovanja iznosi " << min16 << endl
		<< " 17. najslicniji: " << vekidgovornika[prepoznati17] << " pri cemu mera razlikovanja iznosi " << min17 << endl
		<< " 18. najslicniji: " << vekidgovornika[prepoznati18] << " pri cemu mera razlikovanja iznosi " << min18 << endl
		<< " 19. najslicniji: " << vekidgovornika[prepoznati19] << " pri cemu mera razlikovanja iznosi " << min19 << endl
		<< " 20. najslicniji: " << vekidgovornika[prepoznati20] << " pri cemu mera razlikovanja iznosi " << min20 << endl;*/
	/*int broj_stridensa=0,broj_pravilnih=0;
	if(vekidgovornika[prepoznati1]=="Stridens") {broj_stridensa++;} else {broj_pravilnih++;}
	if(vekidgovornika[prepoznati2]=="Stridens") {broj_stridensa++;} else {broj_pravilnih++;}
	if(vekidgovornika[prepoznati3]=="Stridens") {broj_stridensa++;} else {broj_pravilnih++;}
	if(vekidgovornika[prepoznati4]=="Stridens") {broj_stridensa++;} else {broj_pravilnih++;}
	if(vekidgovornika[prepoznati5]=="Stridens") {broj_stridensa++;} else {broj_pravilnih++;}
	if(broj_stridensa>broj_pravilnih){cout<<"U snimku je prisutan stridens"<<endl;} else {cout<<"Izgovor u snimku je pravilan"<<endl;}
	if(vekidgovornika[prepoznati1]=="StevanJokic" | vekidgovornika[prepoznati1]=="IvanJokic")
	{system("\"C:\\Program Files\\Winamp\\winamp.exe\"/NEW Vo_Kosi_Da_Ti_Spijam.mp3");}*/
	//{system("\"C:\\Program Files\\Windows Media Player\\wmplayer.exe\"/NEW Vo_Kosi_Da_Ti_Spijam.mp3");}
		return 0;
}

int ns_prepoznavanje_nad_snimkom(char* snimak, int brojgovornika, vector<string> vekidgovornika,
								 vector<double> vekapsodstkovmat_test)
//(ne!)prikazuje 6 najblizih suseda (govornika) posmatranom test snimku
{
	/*int prepoznati1=0,prepoznati2=0,prepoznati3=0,prepoznati4=0,prepoznati5=0,prepoznati6=0; 
	//pocetne vrednosti rednih brojeva: prepoznatog govornika i drugog najslicnijeg govornika posmatranom govorniku
	double min1=vekapsodstkovmat_test[0],min2=min1,min3=min1,min4=min1,min5=min1,min6=min1;
	for(int i=0;i<brojgovornika;i++)  //radi poredjenja sa modelima svih govornika u govornoj bazi
	{
		if (vekapsodstkovmat_test[i]<min1)
		{
			min6=min5,min5=min4,min4=min3,min3=min2,min2=min1,min1=vekapsodstkovmat_test[i];
			prepoznati6=prepoznati5,prepoznati5=prepoznati4,prepoznati4=prepoznati3,
				prepoznati3=prepoznati2,prepoznati2=prepoznati1,prepoznati1=i;
		}
	}*/
	cout<<"\nZa snimak "<<snimak<<" razlikovanja u odnosu na postojece modele su sledeca:"<<endl 
		//<< " 1. najslicniji: " << vekidgovornika[prepoznati1] << " pri cemu mera razlikovanja iznosi " << min1 << endl
		//<< " 2. najslicniji: " << vekidgovornika[prepoznati2] << " pri cemu mera razlikovanja iznosi " << min2 << endl
		//<< " 3. najslicniji: " << vekidgovornika[prepoznati3] << " pri cemu mera razlikovanja iznosi " << min3 << endl
		//<< " 4. najslicniji: " << vekidgovornika[prepoznati4] << " pri cemu mera razlikovanja iznosi " << min4 << endl
		//<< " 5. najslicniji: " << vekidgovornika[prepoznati5] << " pri cemu mera razlikovanja iznosi " << min5 << endl
		//<< " 6. najslicniji: " << vekidgovornika[prepoznati6] << " pri cemu mera razlikovanja iznosi " << min6 << endl;
	    << " u odnosu na " << vekidgovornika[0] << " mera razlikovanja iznosi " << vekapsodstkovmat_test[0] << endl
		<< " u odnosu na " << vekidgovornika[1] << " mera razlikovanja iznosi " << vekapsodstkovmat_test[1] << endl
		<< " u odnosu na " << vekidgovornika[2] << " mera razlikovanja iznosi " << vekapsodstkovmat_test[2] << endl
		<< " u odnosu na " << vekidgovornika[3] << " mera razlikovanja iznosi " << vekapsodstkovmat_test[3] << endl
		<< " u odnosu na " << vekidgovornika[4] << " mera razlikovanja iznosi " << vekapsodstkovmat_test[4] << endl
		<< " u odnosu na " << vekidgovornika[5] << " mera razlikovanja iznosi " << vekapsodstkovmat_test[5] << endl;

	//if(vekidgovornika[prepoznati]=="StevanJokic" | vekidgovornika[prepoznati]=="IvanJokic")
	//{system("\"C:\\Program Files\\Winamp\\winamp.exe\"/NEW Vo_Kosi_Da_Ti_Spijam.mp3");}
		return 0;
}

int upisivanje_modela(int broj_govornika,vector<string> vektor_govornika,int brojMFCCs,
					  vector<vector<vector<double>>> kovmat1,vector<vector<vector<double>>> kovmat2,
					  vector<vector<vector<double>>> kovmat1_fajlovi)
{
	fstream modeli("modeli.txt", ios::out|ios::binary);
	fstream modeli1("modeli1.txt", ios::out|ios::binary);//za modele svakog fajla pri obuci
	for (int i=0;i<broj_govornika;i++)
	{
		modeli << "\n" << vektor_govornika[i] << "\n";  //ubacivanje imena govornika u modeli.txt fajl iza koga slede elementi kovarijansne matrice, skupa za obuku i test fajla
		for (int l=0; l<brojMFCCs; l++)
		{
			for (int m=0; m<brojMFCCs; m++)
			{
				modeli << "KovMat " << l+1 << ", " << m+1 << "   " << kovmat1[i][l][m] << "   " << kovmat2[i][l][m] << "   " << kovmat1[i][l][m]-kovmat2[i][l][m] << "\n";  //radi uporednog prikaza kovmat dobijene pri obuci i za test fajl
			}
		}
	} //i
	modeli.close();
	for(int i=0;i<broj_govornika;i++) 
	{modeli1<<"\n"<<vektor_govornika[i]<<"\n";
	for(int j=0;j<brojMFCCs;j++) 
	{for(int k=0;k<brojMFCCs;k++)
	{for(int l=0;l<=12;l++) //za 13 fajlova koji se koriste za obuku
	{ if(l==0) {modeli1<<"KovMat"<<j+1<<", "<<k+1<<"  ";}
			modeli1<< kovmat1_fajlovi[i+l][j][k] << "  ";
	     //ispis svih elemenata sa istog mesta iz kov.matrica u jednoj vrsti, uporedo
	if(l==12){modeli1<<"\n";/*za prelazak u novi red*/}}}}}
	modeli1.close();
	return 0;
}

vector<vector<vector<double>>> ucitavanje_modela(int broj_govornika,int brojMFCCs,const char *imefajla)
{
	vector<vector<vector<double>>> kovmat1(broj_govornika, vector<vector<double>> (brojMFCCs, vector<double> (brojMFCCs)));   // trening kovarijansne matrice
	//vector<vector<vector<float>>> kovmat2(broj_govornika, vector<vector<float>> (brojMFCCs, vector<float> (brojMFCCs)));
	ifstream modeli;
	modeli.open(imefajla);
	string linija; string prvi_znak; int vrsta,kolona; double vrednost,vrednost_t; int brojac_govornika=-1;
	while (!modeli.eof())
	{
		string oznaka_vrste="",oznaka_kolone="",oznaka_vrednosti="",oznaka_vrednosti_t="";
		getline(modeli, linija); //nakon sto se procita linija potrebno je ustanoviti kakve informacije ona nosi
		                       //da li je prazna, da li sadrzi oznaku govornika ili sadrzi elemente matrice
		                       //zatim vrednosti upisivati na odgovarajuca mesta u matrici kovmat1 i kovmat2
		prvi_znak=linija[0];
		if(prvi_znak=="s" | prvi_znak=="I") {brojac_govornika++;}
		int brojac_praznina=0; string indikacija;
		string::iterator pocetak_vrste,kraj_vrste,pocetak_kolone,kraj_kolone,pocetak_vrednosti,kraj_vrednosti;
		if(prvi_znak=="K") 
		{
			for(int i=0;i<linija.size();i++)
			{
				indikacija=linija[i]; if(indikacija==" "){brojac_praznina++;}
				if((indikacija==" ") && (brojac_praznina==1)) {pocetak_vrste = linija.begin()+i+1;} 
				if(indikacija==",") {kraj_vrste = linija.begin()+i; oznaka_vrste.append(pocetak_vrste,kraj_vrste);
				vrsta=atoi(oznaka_vrste.c_str());} //izmedju prve praznine i zareza je broj reda u matrici
				
				if((indikacija==" ") && (brojac_praznina==2)) {pocetak_kolone = linija.begin()+i+1;}
				if((indikacija==" ") && (brojac_praznina==3)) {kraj_kolone = linija.begin()+i; oznaka_kolone.append(pocetak_kolone,kraj_kolone); 
				kolona=atoi(oznaka_kolone.c_str());} //izmedju druge praznine i trece praznine je broj kolone u matrici
				
				if((indikacija==" ") && brojac_praznina==5) {pocetak_vrednosti = linija.begin()+i+1;}
				if((indikacija==" ") && brojac_praznina==6) {kraj_vrednosti = linija.begin()+i;
				oznaka_vrednosti.append(pocetak_vrednosti,kraj_vrednosti);
				vrednost=atof(oznaka_vrednosti.c_str());} //izmedju pete i seste praznine je vrednost

				//if((indikacija==" ") && brojac_praznina==8) {pocetak_vrednosti = linija.begin()+i+1;}
				//if((indikacija==" ") && brojac_praznina==9) {kraj_vrednosti = linija.begin()+i;
				//oznaka_vrednosti_t.append(pocetak_vrednosti,kraj_vrednosti);
				//vrednost_t=atof(oznaka_vrednosti_t.c_str());} //izmedju osme i devete praznine je test vrednost
			}
			kovmat1[brojac_govornika][vrsta-1][kolona-1]=vrednost;
			//kovmat2[brojac_govornika][vrsta-1][kolona-1]=vrednost_t;
		}
	}
	return kovmat1;
}

vector<vector<vector<double>>> ucitavanje_modela_t(int broj_govornika,int brojMFCCs,const char *imefajla)
{
	//vector<vector<vector<float>>> kovmat1(broj_govornika, vector<vector<float>> (brojMFCCs, vector<float> (brojMFCCs)));   // trening kovarijansne matrice
	vector<vector<vector<double>>> kovmat2(broj_govornika, vector<vector<double>> (brojMFCCs, vector<double> (brojMFCCs)));
	ifstream modeli;
	modeli.open(imefajla);
	string linija; string prvi_znak; int vrsta,kolona; double vrednost,vrednost_t; int brojac_govornika=-1;
	while (!modeli.eof())
	{
		string oznaka_vrste="",oznaka_kolone="",oznaka_vrednosti="",oznaka_vrednosti_t="";
		getline(modeli, linija); //nakon sto se procita linija potrebno je ustanoviti kakve informacije ona nosi
		                       //da li je prazna, da li je sadrzi oznaku govornika ili sadrzi elemente matrice
		                       //zatim vrednosti upisivati na odgovarajuca mesta u matrici kovmat1 i kovmat2
		prvi_znak=linija[0];
		if(prvi_znak=="s" | prvi_znak=="I") {brojac_govornika++;}
		int brojac_praznina=0; string indikacija;
		string::iterator pocetak_vrste,kraj_vrste,pocetak_kolone,kraj_kolone,pocetak_vrednosti,kraj_vrednosti;
		if(prvi_znak=="K") 
		{
			for(int i=0;i<linija.size();i++)
			{
				indikacija=linija[i]; if(indikacija==" "){brojac_praznina++;}
				if((indikacija==" ") && (brojac_praznina==1)) {pocetak_vrste = linija.begin()+i+1;} 
				if(indikacija==",") {kraj_vrste = linija.begin()+i; oznaka_vrste.append(pocetak_vrste,kraj_vrste);
				vrsta=atoi(oznaka_vrste.c_str());} //izmedju prve praznine i zareza je broj reda u matrici
				
				if((indikacija==" ") && (brojac_praznina==2)) {pocetak_kolone = linija.begin()+i+1;}
				if((indikacija==" ") && (brojac_praznina==3)) {kraj_kolone = linija.begin()+i; oznaka_kolone.append(pocetak_kolone,kraj_kolone); 
				kolona=atoi(oznaka_kolone.c_str());} //izmedju druge praznine i trece praznine je broj kolone u matrici
				
				//if((indikacija==" ") && brojac_praznina==5) {pocetak_vrednosti = linija.begin()+i+1;}
				//if((indikacija==" ") && brojac_praznina==6) {kraj_vrednosti = linija.begin()+i;
				//oznaka_vrednosti.append(pocetak_vrednosti,kraj_vrednosti);
				//vrednost=atof(oznaka_vrednosti.c_str());} //izmedju pete i seste praznine je vrednost

				if((indikacija==" ") && brojac_praznina==8) {pocetak_vrednosti = linija.begin()+i+1;}
				if((indikacija==" ") && brojac_praznina==9) {kraj_vrednosti = linija.begin()+i;
				oznaka_vrednosti_t.append(pocetak_vrednosti,kraj_vrednosti);
				vrednost_t=atof(oznaka_vrednosti_t.c_str());} //izmedju osme i devete praznine je test vrednost
			}
			//kovmat1[brojac_govornika][vrsta-1][kolona-1]=vrednost;
			kovmat2[brojac_govornika][vrsta-1][kolona-1]=vrednost_t;
		}
	}
	return kovmat2;
}

vector<vector<vector<double>>> ucitavanje_modela_fajlovi(int broj_govornika,int brojMFCCs,const char *imefajla)
{
	int broj_fajlova=broj_govornika*13;
	vector<vector<vector<double>>> kovmat1_fajlovi(broj_fajlova, vector<vector<double>> (brojMFCCs, vector<double> (brojMFCCs)));
	ifstream modeli;
	modeli.open(imefajla);
	string linija,prvi_znak; int vrsta,kolona; double vrednost; int brojac_govornika=-1;
	while (!modeli.eof())
	{
		string oznaka_vrste="",oznaka_kolone="",oznaka_vrednosti="";
		getline(modeli, linija); //nakon sto se procita linija potrebno je ustanoviti kakve informacije ona nosi
		                       //da li je prazna, da li sadrzi oznaku govornika ili sadrzi elemente matrice
		                       //zatim vrednosti upisivati na odgovarajuca mesta u matrici kovmat1_fajlovi
		prvi_znak=linija[0];
		if(prvi_znak=="s") {brojac_govornika++;}
		int brojac_praznina=0; string indikacija;
		string::iterator pocetak_vrste,kraj_vrste,pocetak_kolone,kraj_kolone,pocetak_vrednosti,kraj_vrednosti;
		if(prvi_znak=="K") 
		{
			for(int i=0;i<strlen(linija.c_str());i++)
			{
				indikacija=linija[i]; if(indikacija==" "){brojac_praznina++;}
				if((indikacija=="t") && (brojac_praznina==0)) {pocetak_vrste = linija.begin()+i+1;} 
				if(indikacija==",") {kraj_vrste = linija.begin()+i; oznaka_vrste.append(pocetak_vrste,kraj_vrste);
				vrsta=atoi(oznaka_vrste.c_str());} //izmedju slova "t" i zareza je broj reda u matrici
				
				if((indikacija==" ") && (brojac_praznina==1)) {pocetak_kolone = linija.begin()+i+1;}
				if((indikacija==" ") && (brojac_praznina==2)) {kraj_kolone = linija.begin()+i; oznaka_kolone.append(pocetak_kolone,kraj_kolone); 
				kolona=atoi(oznaka_kolone.c_str());} //izmedju druge i trece praznine je broj kolone u matrici
				
				for(int j=1;j<=13;j++)
				{
				if((indikacija==" ") && brojac_praznina==2*j+1) {pocetak_vrednosti = linija.begin()+i+1;}
				if((indikacija==" ") && brojac_praznina==2*j+2) {kraj_vrednosti = linija.begin()+i;
				oznaka_vrednosti.append(pocetak_vrednosti,kraj_vrednosti);
				vrednost=atof(oznaka_vrednosti.c_str()); //izmedju 2*j+1 i 2*j+2 praznine je vrednost
				kovmat1_fajlovi[brojac_govornika+j-1][vrsta-1][kolona-1]=vrednost;
				oznaka_vrednosti="";}//inicijalizacija oznake vrednosti
				}
			}
			
		}
	}
	return kovmat1_fajlovi;
}


vector<vector<double>> matrica_razlikovanja (int broj_govornika,int brojMFCCs,vector<vector<vector<double>>> kovmat1,vector<vector<vector<double>>> kovmat2)
{//izracunava ukupno razlikovanje elementa kov. matrice na nivou cele baze
	fstream razlikovanja("razlikovanja.txt", ios::out|ios::binary);
	vector<vector<double>> kov_mat_raz(brojMFCCs,vector<double>(brojMFCCs));
	for (int i=0;i<broj_govornika;i++)
	{
		for(int j=0;j<brojMFCCs;j++)
		{
			for(int k=0;k<brojMFCCs;k++)
			{
				kov_mat_raz[j][k]=kov_mat_raz[j][k]+abs(kovmat1[i][j][k]-kovmat2[i][j][k]);
			}
		}
	} //i
	double min=kov_mat_raz[0][0];
	for (int i=0;i<brojMFCCs;i++)
	{
		for (int j=0; j<brojMFCCs; j++)
		{
			if(kov_mat_raz[i][j]<min) {min=kov_mat_raz[i][j];}
		}
	}
	razlikovanja << "Najmanje razlikovanje iznosi: " << min << endl;

	for (int j=0; j<brojMFCCs; j++)
	{
		for (int k=0; k<brojMFCCs; k++)
		{
				razlikovanja << "RazKovMat " << j+1 << ", " << k+1 << "   " << kov_mat_raz[j][k] << "   RelRazKovMat " << kov_mat_raz[j][k]/min << "\n";  //prikaza razlikovanja kov_mat po pojedinim elementima
		}
	}
	razlikovanja.close();
	return kov_mat_raz;
}

int pravljenje_wav_fajla (vector<double> odbirci,const char *ime)
{
	fstream snimak(ime, ios::out|ios::binary);
	unsigned char b;  //velicine je jednog bajta, upisuje se bajt po bajt
	                  //http://www.cplusplus.com/forum/beginner/63396/
	b=0x52;snimak<<b;b=0x49;snimak<<b;b=0x46;snimak<<b;b=0x46;snimak<<b; //"RIFF"

	int broj_bajtova=(odbirci.size())*2+44-8;//treba zapisati u 4 bajta u little endian formatu
	b=broj_bajtova;//posto je b velicine jedan bajt prepisuje se samo prvi levi bajt vrednosti broja_bajtova
	snimak<<b;
	b=broj_bajtova>>8;snimak<<b;//zatim se cita drugi bajt, prvo se prebaci u krajnju levu poziciju
	b=broj_bajtova>>16;snimak<<b;b=broj_bajtova>>24;snimak<<b;//treci i cetvrti bajt da bi se u datoteku 
	                             //upisali u poretku od najmanje znacajnog ka nejvise znacajnom bajtu

	b=0x57;snimak<<b;b=0x41;snimak<<b;b=0x56;snimak<<b;b=0x45;snimak<<b; //"WAVE"
	b=0x66;snimak<<b;b=0x6D;snimak<<b;b=0x74;snimak<<b;b=0x20;snimak<<b; //"fmt "

	b=0x10;snimak<<b;b=0x00;snimak<<b;b=0x00;snimak<<b;b=0x00;snimak<<b; //velicina "fmt " = 16 bajta
	b=0x01;snimak<<b;b=0x00;snimak<<b; // PCM
	b=0x01;snimak<<b;b=0x00;snimak<<b; // mono
	b=0x22;snimak<<b;b=0x56;snimak<<b;b=0x00;snimak<<b;b=0x00;snimak<<b; // fs=22050Hz
	b=0x44;snimak<<b;b=0xAC;snimak<<b;b=0x00;snimak<<b;b=0x00;snimak<<b; // bajt brzina
	b=0x02;snimak<<b;b=0x00;snimak<<b; // broj bajta po odbirku
	b=0x10;snimak<<b;b=0x00;snimak<<b; // broj bita po odbirku (16)
	b=0x64;snimak<<b;b=0x61;snimak<<b;b=0x74;snimak<<b;b=0x61;snimak<<b; // "data"
	int broj_bajtova_podataka=(odbirci.size())*2; // ovde mono zapis sa dva bajta po odbirku

	b=broj_bajtova_podataka;snimak<<b;b=broj_bajtova_podataka>>8;snimak<<b;
	b=broj_bajtova_podataka>>16;snimak<<b;b=broj_bajtova_podataka>>24;snimak<<b;

	int vrednost;
	for(int i=0;i<odbirci.size();i++) // upis vrednosti signala
	{
		if(odbirci[i]>=0){vrednost=odbirci[i]*32767.0;b=vrednost;snimak<<b;b=vrednost>>8;snimak<<b;}
		if(odbirci[i]<0){vrednost=odbirci[i]*32768.0+65536;b=vrednost;snimak<<b;b=vrednost>>8;snimak<<b;}
	}

	snimak.close();
	return 0;
}

void SaveWaveFile(char* ime_govornika, int nacin_rada)
{
	MMCKINFO ChunkInfo;
	MMCKINFO FormatChunkInfo;
	MMCKINFO DataChunkInfo;

	//string a=string(ime_govornika),b=".wav",c=a+b; cout <<c<<endl; const char *ime=c.c_str(); char *ime_wav_fajla=const_cast<char *> (ime);
	char *ime_wav_fajla; if(nacin_rada==2){ime_wav_fajla="testiranje.wav";} if(nacin_rada==1){ime_wav_fajla="obucavanje.wav";}
	//if(obuka==3){string a=string(ime_govornika),b="_obuka.wav",c=a+b; cout <<c<<endl; const char *ime=c.c_str(); ime_wav_fajla=const_cast<char *> (ime);}
	//if(obuka==2){string a=string(ime_govornika),b="_test.wav",c=a+b; cout <<c<<endl; const char *ime=c.c_str(); ime_wav_fajla=const_cast<char *> (ime);}

	HMMIO handle = mmioOpen(
		ime_wav_fajla, 0, MMIO_CREATE | MMIO_WRITE);     //"testiranje.wav"
	if (!handle) {
		MessageBox(0, "Error creating file.", "Error Message", 0);
		return;
	}

	memset(&ChunkInfo, 0, sizeof(MMCKINFO));
	ChunkInfo.fccType = mmioStringToFOURCC("WAVE", 0);
	DWORD Res = mmioCreateChunk(
	handle, &ChunkInfo, MMIO_CREATERIFF);
	CheckMMIOError(Res);

	FormatChunkInfo.ckid = mmioStringToFOURCC("fmt ", 0);
	FormatChunkInfo.cksize = sizeof(WAVEFORMATEX);
	Res = mmioCreateChunk(handle, &FormatChunkInfo, 0);
	CheckMMIOError(Res);
	// Write the wave format data.
	mmioWrite(handle, (char*)&pFormat, sizeof(pFormat));

	Res = mmioAscend(handle, &FormatChunkInfo, 0);
	CheckMMIOError(Res);
	DataChunkInfo.ckid = mmioStringToFOURCC("data", 0);
	DWORD DataSize = WaveInHdr.dwBytesRecorded;
	DataChunkInfo.cksize = DataSize;
	Res = mmioCreateChunk(handle, &DataChunkInfo, 0);
	CheckMMIOError(Res);

	mmioWrite(handle, (char*)WaveInHdr.lpData, DataSize);
	// Ascend out of the data chunk.
	mmioAscend(handle, &DataChunkInfo, 0);

	// Ascend out of the RIFF chunk (the main chunk). Failure to do 
	// this will result in a file that is unreadable by Windows95
	// Sound Recorder.
	mmioAscend(handle, &ChunkInfo, 0);
	mmioClose(handle, 0);
}

void CheckMMIOError(DWORD code)
{
	if (code == 0)
	{
		return;
	}
	
	char buff[256];
	wsprintf(buff, "MMIO Error. Error Code: %d", code);
	MessageBox(NULL,buff, "MMIO Error", 0);
}

int main(int argc, char* argv[])   //broj argumenata main funkcije, uracunavajuci i ime samog programa
{  
	int nacin_rada; //obucavanje na govornika (1) ili prepoznavanje govornika u odredjenom snimku (2)
	nacin_rada = atoi(argv[1]); // 1-nad govornom bazom i prepoznavanje nad odredjenim snimkom, 2-pracenje govornika;
	char* snimak_za_obuku;
	snimak_za_obuku = argv[2];
	char* ime_govornika;
	ime_govornika = argv[3];
	char* snimak_za_prepoznavanje;
	snimak_za_prepoznavanje = argv[4];
	char* govornici;
	govornici = argv[5];
	char* prag;
	prag = argv[6];
	//char* br_dodatnih_obelezja;
	//br_dodatnih_obelezja = argv[7];
	char* ind_post_prvog_dodatnog_obelezja; //indikator postojanja prvog dodatnog obelezja, 0 ako se ne uzima, 1 ako se uzima u obzir
	ind_post_prvog_dodatnog_obelezja = argv[7];
	char* ind_post_drugog_dodatnog_obelezja; //indikator postojanja drugog dodatnog obelezja, ovi indikatori su uvedeni da bi se pratilo 
	ind_post_drugog_dodatnog_obelezja = argv[8];// kako i drugo dodatno obelezje ponaosob utice na tacnost prepoznavanja
	char* ind_post_treceg_dodatnog_obelezja;
	ind_post_treceg_dodatnog_obelezja = argv[9];

	char* oblik_kr_opsega="sigm";//"e_d";//"sigm" //definisani su eksponencijalni kriticni osezi
	int brkropsega=22; //20, 22
	//cout<<"Prepoznavac: 20 e_d kriticnih opsega,s=2, 18 MFCCs."<<endl;
	cout<<"\nPrepoznavac prilagodjen radu nad signalom cija je fs=44100 Hz, analizira se prozor od 1024 odbirka (FFT u 1024 tacke)."<<endl;

	//float N = 512.0; //broj tacaka signala za koje se racuna FFT  25ms~512, za fs=22050Hz
	float N = 1024.0;//za fs=44100Hz
	//int pomerajprozora=184;   //8.33ms   184, za fs=22050Hz
	int pomerajprozora=368; //za fs=44100Hz
	int brojMFCCs=21;//18, 22; //definicija broja MFCCs
	//pri racunanju MFCCs pomocu 22 opsega ucestanosti treba za brojMFCCs uzeti 21 posto bi 22. koeficijent bio nula
	//22 koeficijenta bi se mogla ukupno racunati sa nultim ali posto sada ne racunam nulti koeficijent
	//vektor MFCCs ce imati 21 koeficijent
	//cout<<"Prepoznavac: "<<brkropsega<<" e_d kriticnih opsega, s=2, "<<brojMFCCs<< " MFCCs"<<endl;
	cout<<"\nPrepoznavac: "<<brkropsega<<" "<<oblik_kr_opsega<<" kriticna opsega, "<<brojMFCCs<<" MFCCs"<<endl;
	//int broj_dodatnih_obelezja = atoi(br_dodatnih_obelezja);
	int prvo_dodatno_obelezje = atoi(ind_post_prvog_dodatnog_obelezja);
	int drugo_dodatno_obelezje = atoi(ind_post_drugog_dodatnog_obelezja);
	int trece_dodatno_obelezje = atoi(ind_post_treceg_dodatnog_obelezja);

	int broj_dodatnih_obelezja = prvo_dodatno_obelezje + drugo_dodatno_obelezje + trece_dodatno_obelezje;
	
	int broj_obelezja = brojMFCCs + broj_dodatnih_obelezja;
	
	int broj_energetskih_obelezja = 12;
	
	if(prvo_dodatno_obelezje==0 && drugo_dodatno_obelezje==0)   //ispisi o tome koja se obelezja koriste
	{cout<<"\nObelezja: "<<brojMFCCs<<" MFCCs"<<endl;}
	if(prvo_dodatno_obelezje==1 && drugo_dodatno_obelezje==0)
	{cout<<"\nObelezja: "<<brojMFCCs<<" MFCCs + e1, ukupno "<<broj_obelezja<<" obelezja"<<endl;}
	if(prvo_dodatno_obelezje==0 && drugo_dodatno_obelezje==1)
	{cout<<"\nObelezja: "<<brojMFCCs<<" MFCCs + e2, ukupno "<<broj_obelezja<<" obelezja"<<endl;}
	if(prvo_dodatno_obelezje==1 && drugo_dodatno_obelezje==1)
	{cout<<"\nObelezja: "<<brojMFCCs<<" MFCCs + e1 + e2, ukupno "<<broj_obelezja<<" obelezja"<<endl;}
	if(prvo_dodatno_obelezje==1 && drugo_dodatno_obelezje==1 && trece_dodatno_obelezje==1)
	{cout<<"\nObelezja: "<<brojMFCCs<<" MFCCs + e1 + e2 + e3, ukupno "<<broj_obelezja<<" obelezja"<<endl;}

	if(nacin_rada==3) {
		cout << "Vektor obelezja ukupno ima " << broj_obelezja << " obelezja" << endl; }
	
	if(nacin_rada==7 || nacin_rada==8) {
		cout<<"Vektor obelezja ima "<<broj_energetskih_obelezja<<" obelezja"<<endl; }
	
	vector<float> prozor(N); //nesto manje od broja odbiraka u prozoru duzine 25ms, fs=22050Hz, radi uskladjivanja sa duzinom potrebnom za FFT (2^n=512, tacno bi bilo 551 tacka)
	for (int j=0; j<N; j++)
	{
		//prozor[j] = 0.54-0.46*cos(2*PI*j/(N-1));  //standardan Hamingov prozor
		prozor[j] = 0.5*(1-cos(2*PI*j/(N-1))); //Hanov prozor
		//prozor[j] = 0.42659 - 0.49656*cos(2*PI*j/(N-1)) + 0.076849*cos(4*PI*j/(N-1)); //Blekmanov prozor
		//prozor[j] = 3.20-3.10*cos(2*PI*j/(N-1));  //eksperimentalni prozor
	}

	vector<int> obrnutired(N);
	for (int i=0; i<N; i++)
	{
		//obrnutired[i] = ((i>>8) | ((i>>6)&2) | ((i>>4)&4) | ((i>>2)&8) | (i&16) | ((i<<2)&32) | ((i<<4)&64) | ((i<<6)&128) | ((i<<8)&256));   // N=512
		obrnutired[i] = ((i>>9) | ((i>>7)&2) | ((i>>5)&4) | ((i>>3)&8) | ((i>>1)&16) | ((i<<1)&32) | ((i<<3)&64) | ((i<<5)&128) | ((i<<7)&256) | ((i<<9)&512)); // N=1024
	}

	vector<float> costab (N/2);
	vector<float> sintab (N/2);
	for (int i=0; i<(N/2); i++)
	{
		costab[i] = cos(2*PI*i/N);
		sintab[i] = sin(2*PI*i/N);  //racunanje vektora rotacionih faktora za FFT u N tacaka
	}

	vector<vector<float>> costabMFCC (brojMFCCs, vector<float> (brkropsega));  //tabela (matrica) vrednosti kosinusa koje se koriste u DCT za racunanje MFCC
	for (int j=0; j<brojMFCCs; j++)  
	{
		for (int i=0; i<brkropsega; i++)  
		{
			float argument;
			//argument = j*(i+1.0-0.5);  //ukoliko se u obzir uzima nulti i prvih brojMFCCs-1 MFCCs
			argument=(j+1)*(i+1.0-0.5)*PI/22.0;//*PI/22.0 //ukoliko se ne uzima u obzir nulti MFCC, nego prvih brojMFCC MFCCs, i+1 stoji zato sto u jednakosti za racunanje MFCCs i pocinje od 1
			costabMFCC[j][i] = cos(argument);
		}
	}
	
		
	if(nacin_rada==1 | nacin_rada==2) //pravljenje wav snimka sa mikrofona za obuku ili test - pocetak
	{
		//Declare local varibles
	int samplesperchan = 16; 
	int sampleRate = 22050;//44100;
	int *waveIn;

/*	cout << "*********************************************\n";
	cout << "Configuring the Sound Hardware:\n";
	cout << "*********************************************\n";
	cout << "Enter the number of Samples/Channel:\n";
	cin >> samplesperchan;
	cout << "Enter the Sampling Rate:\n";
	cin >> sampleRate; */

	pFormat.wFormatTag = WAVE_FORMAT_PCM; // simple, uncompressed format
	pFormat.nChannels = 1; // 1=mono, 2=stereo
	pFormat.nSamplesPerSec = sampleRate; // 44100
	pFormat.wBitsPerSample = 16; // 16 for high quality, 8 for telephone-grade
	pFormat.nBlockAlign = pFormat.nChannels*pFormat.wBitsPerSample/8; //broj bajta po odbirku
	//pFormat.nAvgBytesPerSec = pFormat.nChannels*pFormat.wBitsPerSample/8;
	pFormat.nAvgBytesPerSec = sampleRate*pFormat.nChannels*pFormat.wBitsPerSample/8; //broj bajta u sekundi (bajtska brzina)
	pFormat.cbSize=0;

	UINT waveInNumDevs = waveInGetNumDevs();

	if (!waveInNumDevs)
	{
		MessageBoxW(NULL, L"No recording devices found in the system.", NULL, NULL);
		//return;
	}

	result = waveInOpen(&hWaveIn, WAVE_MAPPER, &pFormat, 0L, 0L, WAVE_FORMAT_DIRECT);
	if (result)
	{
		char fault[256];
		waveInGetErrorText(result, fault, 256);
		cout << " result:" << result << endl << fault << endl << WAVE_MAPPER << endl << waveInNumDevs << endl;
		MessageBoxW(NULL,L"Failed to open waveform input device. ",NULL,NULL);		
		//return;
	}
	
	int nSec = NSEC;
	waveIn = new int[sampleRate*nSec];
	WaveInHdr.lpData = (LPSTR)waveIn;
	WaveInHdr.dwBufferLength = sampleRate*nSec*pFormat.nBlockAlign;
	WaveInHdr.dwBytesRecorded=0;
	WaveInHdr.dwUser = 0L;
	WaveInHdr.dwFlags = 0L;
	WaveInHdr.dwLoops = 0L;
	waveInPrepareHeader(hWaveIn, &WaveInHdr, sizeof(WAVEHDR));

	result = waveInAddBuffer(hWaveIn, &WaveInHdr, sizeof(WAVEHDR));
	if (result)
	{
		MessageBoxW(NULL,L"Failed to read block from device",NULL,NULL);
		//return;
	}

	result = waveInStart(hWaveIn);
	if (result)
	{
		MessageBoxW(NULL,L"Failed to start recording",NULL,NULL);
		//return;
	}
	cout << "Snimanje pocinje i traje 10 sekundi...........\n";
	do 
	{
	}while (waveInUnprepareHeader(hWaveIn, &WaveInHdr, sizeof(WAVEHDR))== WAVERR_STILLPLAYING);

	SaveWaveFile(ime_govornika,nacin_rada);

	waveInStop(hWaveIn);

	cout << "Snimanje je zavrseno.\n";

	waveInClose(hWaveIn);

	if (!waveIn) 
	{
		delete[] waveIn;
	}
	}  //pravljenje wav snimka sa mikrofona za obuku ili test - kraj

	int broj_govornika; vector<string> vektor_govornika(12000);//proveriti da li je potreban ovaj vektor posto ga u nekim nacinima
//rada ne koristim npr. za nacin rada 11 a moze da izazove zastoj rada programa ako je dimenzija vektora manja od broja modela 
	//govornika u bazi govornika 
	if(_access("govornici.txt",0)!=-1) //ako nije prazna govorna baza vrsi se prebrojavanje i upis govornika u vektor govornika
	{
		ifstream baza; baza.open(govornici);
		broj_govornika=0; //u sustini ovo ce biti broj modela govornika
		string govornik;
		while (!baza.eof())
		{
			getline(baza, govornik);
			broj_govornika++;
		}
		cout << "\nGovorna baza ima " << broj_govornika << " modela govornika." << endl;
		
		ifstream baza1; baza1.open(govornici);
		for (int i=0; i<broj_govornika; i++)
		{
			getline(baza1, vektor_govornika[i]);  //prepisuje nazive govornika iz govornici.txt u vektor_govornika
		}
	}

	if(nacin_rada==1 | nacin_rada==3) //obucavanje na govornika
	{
		ofstream baza_g;
		if(_access("govornici.txt",0)==-1){baza_g.open("govornici.txt");baza_g<<ime_govornika;baza_g.close();broj_govornika=1;} //ako nema baze govornika i program je prvi put pokrenut
		else{baza_g.open("govornici.txt",std::ios::app);baza_g<<"\n"<<ime_govornika;baza_g.close();broj_govornika++;} //ubacuje ime govornika za koga ce se vrsiti obuka u bazu imena govornika
			
		vector<double> odbirci; vector<double> vek_sr_vr; vector<vector<double>> kov_mat;
		
		odbirci.swap(ucitavanjewavfajla(snimak_za_obuku));
		brojprozora = ((odbirci.size())-N)/pomerajprozora; //racunanje modela snimljenog fajla za obuku
		vector<vector<double>> mfccs(brojprozora, vector<double> (brojMFCCs));
		mfccs.swap(mfcc(N, pomerajprozora, odbirci, oblik_kr_opsega, brkropsega, brojMFCCs, obrnutired, prozor, costab, sintab, costabMFCC));
		vector<vector<double>> vektori_dodatnih_obelezja(brojprozora,vector<double> (3));
		vektori_dodatnih_obelezja.swap(dodatna_obelezja(N, pomerajprozora, odbirci, obrnutired, prozor, costab, sintab,prvo_dodatno_obelezje,drugo_dodatno_obelezje,broj_dodatnih_obelezja));	
		vector<vector<double>> vektori_obelezja(brojprozora, vector<double> (broj_obelezja));
		if(broj_dodatnih_obelezja == 0) {vektori_obelezja.swap(mfccs);}
		if(broj_dodatnih_obelezja != 0) 
		{
			for (int i=0; i<brojprozora; i++)
			{
				for(int j=0; j<broj_obelezja; j++)
				{
					if(j<brojMFCCs) {vektori_obelezja[i][j] = mfccs[i][j];}
					//if(j>=brojMFCCs) {vektori_obelezja[i][j] = vektor_dodatnih_obelezja[i];}
					
					//if(j<brojMFCCs-1) {vektori_obelezja[i][j] = mfccs[i][j];}
					//if(j==brojMFCCs-1) {vektori_obelezja[i][j] = vektor_dodatnih_obelezja[i];}
					
					if(j>=brojMFCCs) {
						if(broj_dodatnih_obelezja==1 || broj_dodatnih_obelezja==2 || broj_dodatnih_obelezja==3)
						{vektori_obelezja[i][brojMFCCs]=vektori_dodatnih_obelezja[i][0];}
						if(broj_dodatnih_obelezja==2 || broj_dodatnih_obelezja==3)
						{vektori_obelezja[i][brojMFCCs+1]=vektori_dodatnih_obelezja[i][1];}
						if(broj_dodatnih_obelezja==3)
						{vektori_obelezja[i][brojMFCCs+2]=vektori_dodatnih_obelezja[i][2];}
					}
				}
			}
		}
		vek_sr_vr.swap(veksrvr(broj_obelezja,vektori_obelezja));
		kov_mat.swap(kovmat(broj_obelezja,vektori_obelezja,vek_sr_vr));
		
		ofstream modeli; modeli.open("modeli.txt", std::ios::app); // model se dodaje na kraj txt baze modela
		modeli << "\n" << "s" << ime_govornika << "\n";  //ubacivanje imena govornika u modeli.txt fajl iza koga slede elementi kovarijansne matrice, skupa za obuku i test fajla
		for (int i=0; i<broj_obelezja; i++)
		{
			for (int j=0; j<broj_obelezja; j++)
			{
				modeli << "KovMat " << i+1 << ", " << j+1 << "   " << kov_mat[i][j] << "   " << kov_mat[i][j] << "   " << kov_mat[i][j]-kov_mat[i][j] << "\n";  //radi uporednog prikaza kovmat dobijene pri obuci i za test fajl
			}
		}
		modeli.close();
		cout<<"\nObuka prepoznavaca na govornika '"<<ime_govornika<<"' je uspesno izvrsena i sada u bazi modela ima "<<broj_govornika<< " modela"<<endl;
	}
	if(nacin_rada==2 | nacin_rada==4) //prepoznavanje nad odredjenim snimkom iz govorne baze 
	{
		if(_access("govornici.txt",0)==-1){cout<<"Baza govornika je prazna, prvo pokrenuti obuku, opcija 1" << endl;}
		else{
		vector<vector<vector<double>>> kovmat1;
		kovmat1.swap(ucitavanje_modela(broj_govornika,broj_obelezja,"modeli.txt"));//umesto brojMFCCs stavljeno opstije broj_obelezja
		
		vector<double> odbirci_test; vector<double> veksrvr_test; vector<vector<double>> kovmat_test;
		
		odbirci_test.swap(ucitavanjewavfajla(snimak_za_prepoznavanje));
		brojprozora = ((odbirci_test.size())-N)/pomerajprozora;
		vector<vector<double>> mfccs_test(brojprozora, vector<double> (brojMFCCs));
		mfccs_test.swap(mfcc(N, pomerajprozora, odbirci_test, oblik_kr_opsega, brkropsega, brojMFCCs, obrnutired, prozor, costab, sintab, costabMFCC));
		vector<vector<double>> vektori_dodatnih_obelezja(brojprozora,vector<double> (3));
		vektori_dodatnih_obelezja.swap(dodatna_obelezja(N, pomerajprozora, odbirci_test, obrnutired, prozor, costab, sintab,prvo_dodatno_obelezje,drugo_dodatno_obelezje,broj_dodatnih_obelezja));	
		vector<vector<double>> vektori_obelezja(brojprozora, vector<double> (broj_obelezja));
		if(broj_dodatnih_obelezja == 0) {vektori_obelezja.swap(mfccs_test);}
		if(broj_dodatnih_obelezja != 0) 
		{
			for (int i=0; i<brojprozora; i++)
			{
				for(int j=0; j<broj_obelezja; j++)
				{
					if(j<brojMFCCs) {vektori_obelezja[i][j] = mfccs_test[i][j];}
					//if(j>=brojMFCCs) {vektori_obelezja[i][j] = vektor_dodatnih_obelezja[i];}
					
					//if(j<brojMFCCs-1) {vektori_obelezja[i][j] = mfccs[i][j];}
					//if(j==brojMFCCs-1) {vektori_obelezja[i][j] = vektor_dodatnih_obelezja[i];}
					
					if(j>=brojMFCCs) {
						if(broj_dodatnih_obelezja==1 || broj_dodatnih_obelezja==2 || broj_dodatnih_obelezja==3)
						{vektori_obelezja[i][brojMFCCs]=vektori_dodatnih_obelezja[i][0];}
						if(broj_dodatnih_obelezja==2 || broj_dodatnih_obelezja==3)
						{vektori_obelezja[i][brojMFCCs+1]=vektori_dodatnih_obelezja[i][1];}
						if(broj_dodatnih_obelezja==3)
						{vektori_obelezja[i][brojMFCCs+2]=vektori_dodatnih_obelezja[i][2];}
					}
				}
			}
		}
		veksrvr_test.swap(veksrvr(broj_obelezja,vektori_obelezja));//umesto brojMFCCs broj_obelezja,vektori_obelezja umesto mfccs_test
		kovmat_test.swap(kovmat(broj_obelezja,vektori_obelezja,veksrvr_test));//umesto brojMFCCs broj_obelezja; vektori_obelezja

		vector<double> vekapsodstkovmat_test(broj_govornika);   //vektor razlikovanja modela test govora i modela dobijenih prilikom obuke
	
		for(int i=0; i<broj_govornika;i++)
		{
			vekapsodstkovmat_test[i] = (apsodstupanjekovmat(broj_obelezja, kovmat1[i], kovmat_test));//broj_obelezja umesto brojMFCCs
		}

		prepoznavanje_nad_snimkom(snimak_za_prepoznavanje, broj_govornika, vektor_govornika, vekapsodstkovmat_test);
		//knn_prepoznavanje_nad_snimkom(snimak_za_prepoznavanje, broj_govornika, vektor_govornika, vekapsodstkovmat_test);
		//ns_prepoznavanje_nad_snimkom(snimak_za_prepoznavanje, broj_govornika, vektor_govornika, vekapsodstkovmat_test);
/*
		//mogla bi se ovde ubaciti petlja po broju prozora koja bi za svaki ram test snimka racunala vektor kmdfft i ispisivala maksimum u spektru
		vector<double> odbirci1(1024); vector<double> spektar;
		for(int i=0;i<brojprozora;i++)
		{
			for(int j=0;j<1024;j++)//ovde sada ide izracunavanje po okvirima signala
			{
				odbirci1[j]=odbirci_test[i*pomerajprozora+j];//popunjavanje vektora odbiraka
			}
			spektar.swap(kmfft(N,odbirci1,obrnutired,prozor,costab,sintab));
			cout << "\n Prozor " << i << endl;
			spektralni_maksimum(N,spektar);//odredjivanje maksimuma u spektru
		}//ovo bi se moglo iskoristiti za izracunavanje novog dodatnog obelezja
*/		
		}//kraj za else
	} //nacin_rada=2

	if(nacin_rada==5)//deo koji vrsi Proveravanje govornika
		//prepoznavac treba da na osnovu liste govornika u bazi, za svakog trenutno posmatranog govornika svi ostali se mogu smatrati uljezima
		//za svaki referentni model koji odgovara najduzem snimku razmatranog govornika i postavljen prag odluke
		//na osnovu ostalih snimaka preostalih govornika iz govorne baze odredi verovatnocu pogresnog prihvatanja i da
		//na osnovu preostalih snimaka za razmatranog govornika odredi verovatnocu pogresnog odbacivanja
		//ovaj postupak treba da uradi za svih 20 snimaka od 20 izabranih govornika i da
		//kao rezultat ispise kolike su verovatnoce pogresnog prihvatanja i odbacivanja za postavljenu vrednost praga
		//Dakle ulazni parametri ove funkcije bi bili:
		//vrednost praga odlucivanja,
		//lista referentnih govornika, moze biti i kao fajl govornici.txt.
		//Prepoznavac ucitava referentni fajl za odgovarajuceg iz liste od 20 govornika, formira ref. model
		//formira modele za sve ostale snimke u govornoj bazi, i svoje i drugih govornika iz baze, i prati razliku
		//U sustini potrebno je isprogramirati da prepoznavac moze na osnovu imena govornika da zna koji ce fajl uzeti za obuku
		//i koje ce fajlove posmatrati pri testiranju,
		//u sustini to se svodi na analizu niza znakova u nazivu svakog fajla na koji naidje
	{
		double prag_odluke = atof(prag);
	//za pogresno odbacivanje kao test fajlovi se posmatraju oni fajlovi koji imaju u svom imenu oznaku referentnog govornika
	//znaci u prvom koraku se cita i-ta linija sa imenom govornika
	//i onda se prati razlika u modelima test fajlova koji imaju taj id u svom nazivu, u odnosu na referentni model
	//znaci treba da postoji i listing svih fajlova u govornoj bazi
	int broj_govornika=36; //broj poznatih, referentnih, govornika
	int broj_snimaka_po_govorniku=37;
	int ukupan_broj_snimaka = broj_govornika*broj_snimaka_po_govorniku; //36*37=1332
	vector<string> fajlovi(ukupan_broj_snimaka); //vektor stringova u koji ce se upisati nazivi svih snimaka u bazi govornika
	ifstream sadrzaj_baze;  sadrzaj_baze.open("Fajlovi.txt"); //pravljenje spiska wav fajlova: dir /b /s *.wav > list.txt
	for(int i=0; i<ukupan_broj_snimaka; i++)
	{
		getline(sadrzaj_baze, fajlovi[i]); //prebacuje nazive iz Fajlovi.txt u vektor fajlova
	}
	vector<string> poznati_govornici(broj_govornika);//vektor stringova za poznate govornike
	ifstream govornici;  govornici.open("PoznatiGovornici.txt");
	for(int i=0;i<(broj_govornika);i++)
	{getline(govornici, poznati_govornici[i]);}//prebacuju se imena poznatih govornika iz .txt fajla preko toka (strima) govornici u vektor govornika
	
	int broj_pogresnih_odbacivanja = 0; int broj_pogresnih_prihvatanja = 0;
	double p_pogresnog_odbacivanja=0.0; double p_pogresnog_prihvatanja=0.0;
	int br_testova_pogr_odb = 0; int br_testova_pogr_prihv = 0; int redni_broj_testa=0;

	int brso = 9;// 1 ili 9; brso - broj snimaka koji se koriste za obuku
	vector<string> nastavci(brso); nastavci[0]="_s15_solo.wav";nastavci[1]="_s33_solo.wav";nastavci[2]="_s01_solo.wav";nastavci[3]="_s07_solo.wav";
	nastavci[4]="_s05_solo.wav";nastavci[5]="_s31_solo.wav";nastavci[6]="_s18_solo.wav";nastavci[7]="_s20_solo.wav";nastavci[8]="_s25_solo.wav";
	vector<string> referentni_snimci(brso);//za cuvanje naziva referentnih snimaka
	
	for(int i=0;i<broj_govornika;i++) {//deo programa u kom se vrse poredjenja, za svakog govornika ponaosob
		cout<<"\n\n\nTestitanje za govornika "<<poznati_govornici[i]<<endl;
		
		string referentni_govornik = poznati_govornici[i];
		
		//string referentni_snimak = referentni_govornik + "_s15_solo.wav";//_s08_solo.wav
		//cout << "Referentni snimak " << referentni_snimak << endl;
		
		   //vector<vector<vector<double>>> kovmat1(broj_govornika, vector<vector<double>> (brojMFCCs, vector<double> (brojMFCCs)));
		
		vector<vector<vector<double>>> vek_kovmat (brso,vector<vector<double>>(broj_obelezja,vector<double>(broj_obelezja)));//vektor matrica modela za svakog govornika se iznova popunjava
		
		for(int rbrn=0;rbrn<brso;rbrn++) //rbrn - redni broj nastavka 
		{
			referentni_snimci[rbrn] = referentni_govornik + nastavci[rbrn]; 

	vector<double> odbirci;
	odbirci.swap(ucitavanjewavfajla(referentni_snimci[rbrn].c_str())); //.c_str() za konverziju stringa u char* ali ne moze u const char
	brojprozora = ((odbirci.size())-N)/pomerajprozora;
	vector<vector<double>> mfccs_ref(brojprozora, vector<double> (brojMFCCs));
	mfccs_ref.swap(mfcc(N, pomerajprozora, odbirci, oblik_kr_opsega, brkropsega, brojMFCCs, obrnutired, prozor, costab, sintab, costabMFCC));
	
	vector<vector<double>> vektori_dodatnih_obelezja_ref(brojprozora,vector<double> (3));
		vektori_dodatnih_obelezja_ref.swap(dodatna_obelezja(N, pomerajprozora, odbirci, obrnutired, prozor, costab, sintab,prvo_dodatno_obelezje,drugo_dodatno_obelezje,broj_dodatnih_obelezja));
		vector<vector<double>> vektori_obelezja_ref(brojprozora, vector<double> (broj_obelezja));
		if(broj_obelezja == brojMFCCs) {vektori_obelezja_ref.swap(mfccs_ref);}
		if(broj_obelezja > brojMFCCs) 
		{
			for (int i=0; i<brojprozora; i++)
			{
				for(int j=0; j<broj_obelezja; j++)
				{
					if(j<brojMFCCs) {vektori_obelezja_ref[i][j] = mfccs_ref[i][j];}
					//if(j>=brojMFCCs) {vektori_obelezja[i][j] = vektor_dodatnih_obelezja[i];}
					
					//if(j<brojMFCCs-1) {vektori_obelezja[i][j] = mfccs_test[i][j];}
					//if(j==brojMFCCs-1) {vektori_obelezja[i][j] = vektor_dodatnih_obelezja[i];}
					
					if(j>=brojMFCCs) {
						if(broj_dodatnih_obelezja==1 || broj_dodatnih_obelezja==2 || broj_dodatnih_obelezja==3)
						{vektori_obelezja_ref[i][brojMFCCs]=vektori_dodatnih_obelezja_ref[i][0];}
						if(broj_dodatnih_obelezja==2 || broj_dodatnih_obelezja==3)
						{vektori_obelezja_ref[i][brojMFCCs+1]=vektori_dodatnih_obelezja_ref[i][1];}
						if(broj_dodatnih_obelezja==3)
						{vektori_obelezja_ref[i][brojMFCCs+2]=vektori_dodatnih_obelezja_ref[i][2];}
					}
				}
			}
		}

	vector<double> veksrvr_ref(broj_obelezja);
	veksrvr_ref.swap(veksrvr(broj_obelezja,vektori_obelezja_ref));
	vector<vector<double>> kovmat_ref (broj_obelezja,vector<double> (broj_obelezja));
	kovmat_ref.swap(kovmat(broj_obelezja,vektori_obelezja_ref,veksrvr_ref));  //kovmat zvucnog segmenta za obuku

	vek_kovmat[rbrn].swap(kovmat_ref);//odgovarajuci referentni model se upisuje u matricu modela
		}//rbrn, zavrsen upis svih referentnih matrica, modela, za razmatranog govornika
		
	for(int j=0;j<1332;j++) {//sada se prolazi kroz sve snimke u govornoj bazi i na osnovu govornika kom pripadaju radi se test pogresnog odbacivanja ili pogresnog prihvatanja
		string test_snimak = fajlovi[j];
		//iz testirenja ce biti izuzeti i dugacki snimci sa oznakama f01, f02, f03 i f04
		string oznaka; oznaka = "";
		string::iterator pocetak1 = test_snimak.begin()+6;   // http://anaturb.net/C/string_exapm.htm
		string::iterator kraj1 = test_snimak.begin()+9;  //ove oznake se nalaze na 7., 8. i 9. mestu u nazivu snimka
		oznaka.append(pocetak1,kraj1);

		//if( (test_snimak != referentni_snimci[0]) && (oznaka != "f01") && (oznaka != "f02") && (oznaka != "f03") && (oznaka != "f04") ) //ako se koristi samo jedan referentni snimak po govorniku
	    if( (test_snimak != referentni_snimci[0]) && (test_snimak != referentni_snimci[1]) && (test_snimak != referentni_snimci[2]) &&
			(test_snimak != referentni_snimci[3]) && (test_snimak != referentni_snimci[4]) && (test_snimak != referentni_snimci[5]) &&
			(test_snimak != referentni_snimci[6]) && (test_snimak != referentni_snimci[7]) && (test_snimak != referentni_snimci[8]) &&
			(oznaka != "f01") && (oznaka != "f02") && (oznaka != "f03") && (oznaka != "f04") )
		{//da se ne bi vrsilo poredjenje sa istim snimkom koji se koristio i pri obuci i sa dugackim snimcima
		
			vector<double> odbirci_test;
			odbirci_test.swap(ucitavanjewavfajla(test_snimak.c_str()));
			brojprozora = ((odbirci_test.size())-N)/pomerajprozora;
			vector<vector<double>> mfccs_test(brojprozora, vector<double> (brojMFCCs));
			mfccs_test.swap(mfcc(N, pomerajprozora, odbirci_test, oblik_kr_opsega, brkropsega, brojMFCCs, obrnutired, prozor, costab, sintab, costabMFCC));
			
			vector<vector<double>> vektori_dodatnih_obelezja_test(brojprozora,vector<double> (3));
		vektori_dodatnih_obelezja_test.swap(dodatna_obelezja(N, pomerajprozora, odbirci_test, obrnutired, prozor, costab, sintab,prvo_dodatno_obelezje,drugo_dodatno_obelezje,broj_dodatnih_obelezja));
		vector<vector<double>> vektori_obelezja_test(brojprozora, vector<double> (broj_obelezja));
		if(broj_obelezja == brojMFCCs) {vektori_obelezja_test.swap(mfccs_test);}
		if(broj_obelezja > brojMFCCs) 
		{
			for (int i=0; i<brojprozora; i++)
			{
				for(int j=0; j<broj_obelezja; j++)
				{
					if(j<brojMFCCs) {vektori_obelezja_test[i][j] = mfccs_test[i][j];}
					//if(j>=brojMFCCs) {vektori_obelezja[i][j] = vektor_dodatnih_obelezja[i];}
					
					//if(j<brojMFCCs-1) {vektori_obelezja[i][j] = mfccs_test[i][j];}
					//if(j==brojMFCCs-1) {vektori_obelezja[i][j] = vektor_dodatnih_obelezja[i];}
					
					if(j>=brojMFCCs) {
						if(broj_dodatnih_obelezja==1 || broj_dodatnih_obelezja==2 || broj_dodatnih_obelezja==3)
						{vektori_obelezja_test[i][brojMFCCs]=vektori_dodatnih_obelezja_test[i][0];}
						if(broj_dodatnih_obelezja==2 || broj_dodatnih_obelezja==3)
						{vektori_obelezja_test[i][brojMFCCs+1]=vektori_dodatnih_obelezja_test[i][1];}
						if(broj_dodatnih_obelezja==3)
						{vektori_obelezja_test[i][brojMFCCs+2]=vektori_dodatnih_obelezja_test[i][2];}
					}
				}
			}
		}
			
			vector<double> veksrvr_test(broj_obelezja);
			veksrvr_test.swap(veksrvr(broj_obelezja,vektori_obelezja_test));
			vector<vector<double>> kovmat_test (broj_obelezja,vector<double> (broj_obelezja));
			kovmat_test.swap(kovmat(broj_obelezja,vektori_obelezja_test,veksrvr_test));  //kovmat test zvucnog segmenta

	//poredjenje modela
	//double mera_razlikovanja = 0.0;  redni_broj_testa++;
	vector<double> mere_razlikovanja(brso); redni_broj_testa++;//vektor razlikovanja u odnosu na svaki od brso (9) modela za razmatranog govornika

	//mera_razlikovanja = (apsodstupanjekovmat(broj_obelezja, kovmat_s, kovmat_test));
	for(int r=0;r<brso;r++) { mere_razlikovanja[r] = apsodstupanjekovmat(broj_obelezja, vek_kovmat[r], kovmat_test); }
	
	int duzina_imena = referentni_govornik.size(); //duzina imena referentnog govornika
			
			string poredbeni; poredbeni = ""; //radi utvrdjivanja na osnovu test snimka da li fajl pripada ili ne pripada razmatranom poznatom govorniku
			//te da li se stoga odredjuje verovatnoca pogresnog odbacivanja ili prihvatanja
			string::iterator pocetak = test_snimak.begin();   // http://anaturb.net/C/string_exapm.htm
			string::iterator kraj = test_snimak.begin()+duzina_imena;
			poredbeni.append(pocetak,kraj);//iz naziva test snimka se vadi string s pocetka naziva koji sadrzi onoliko znakova kolika je duzina_imena
              //na osnovu imena govornika kome pripada test snimak, nalazi se u stringu poredbeni, i na osnovu imena posmatranog referentnog govornika
			//odredjuje se da li se posmatra pogresno odbacivanje ukoliko su ove dve sadrzine iste ili pogresno prihvatanje ukoliko su razlicite
			
			if( (poredbeni == referentni_govornik) && (test_snimak[duzina_imena]=='_')) {//odredjuje se broj pogresnih odbacivanja
				br_testova_pogr_odb++;
				
				//if(mera_razlikovanja>prag_odluke) {broj_pogresnih_odbacivanja+=1; }
				
				int kontrolna_suma=0;  //nacin da se izvrse sva poredjenja 
				for(int indikator=0;indikator<brso;indikator++) { //pracenjem ukupne vrednosti kontrolne sume
					if(mere_razlikovanja[indikator]>prag_odluke) {kontrolna_suma+=1;} } //nacin da se obezbedi logicko i za sva poredjenja
				if(kontrolna_suma==brso) {broj_pogresnih_odbacivanja+=1;}//ako su sva razlikovanja bila veca od praga tada treba povecati broj pogresnih odbacivanja za 1
				
				/*if ( (mere_razlikovanja[0]>prag_odluke) && (mere_razlikovanja[1]>prag_odluke) && (mere_razlikovanja[2]>prag_odluke)
					&& (mere_razlikovanja[3]>prag_odluke) && (mere_razlikovanja[4]>prag_odluke) && (mere_razlikovanja[5]>prag_odluke)
					&& (mere_razlikovanja[6]>prag_odluke) && (mere_razlikovanja[7]>prag_odluke) && (mere_razlikovanja[8]>prag_odluke) )
				{broj_pogresnih_odbacivanja+=1;}*/
				
					//cout<<"Test "<<redni_broj_testa<<" /42732 "<<"  Referentni snimak "<<referentni_snimak<<"  Test snimak "<<test_snimak<<"  Mera razlikovanja "<<mera_razlikovanja<<"  Broj pogresnih odbacivanja "<<broj_pogresnih_odbacivanja<<" / "<<br_testova_pogr_odb<<" /(1152)"<<endl;}
			cout<<"Test "<<redni_broj_testa<<" /42444 za referentnog govornika "<< referentni_govornik <<"  Test snimak "<<test_snimak<<
				"  Broj pogresnih odbacivanja "<<broj_pogresnih_odbacivanja<<" / "<<br_testova_pogr_odb<<" /(864)"<<endl;}
				//cout<<"Test snimak "<<test_snimak<<" Poredbeni "<<poredbeni<<" Mera razlikovanja "<<mera_razlikovanja<<"   Broj pogresnih odbacivanja "<<broj_pogresnih_odbacivanja<<endl;}   //uvodi se dodatna provera da li je sledeci znak u imenu _ 
                                                                                  //da bi se u obzir uzimali samo fajlovi od posmatranog govornika
			if(poredbeni != referentni_govornik) {//odredjuje se broj pogresnih prihvatanja
				br_testova_pogr_prihv++;
				
				//if(mera_razlikovanja<prag_odluke) {broj_pogresnih_prihvatanja+=1; }

				int kontrolna_suma1=0;  //nacin da se izvrse sva poredjenja 
				for(int indikator1=0;indikator1<brso;indikator1++) { //pracenjem ukupne vrednosti kontrolne sume
					if(mere_razlikovanja[indikator1]<prag_odluke) {kontrolna_suma1+=1;} } //nacin da se obezbedi logicko ili za sva poredjenja
				if(kontrolna_suma1>0) {broj_pogresnih_prihvatanja+=1;}//ako je bar jedno razlikovanje vece od praga tada treba povecati broj pogresnih prihvatanja za 1

				/*if ( (mere_razlikovanja[0]<prag_odluke) || (mere_razlikovanja[1]<prag_odluke) || (mere_razlikovanja[2]<prag_odluke)
					|| (mere_razlikovanja[3]<prag_odluke) || (mere_razlikovanja[4]<prag_odluke) || (mere_razlikovanja[5]<prag_odluke)
					|| (mere_razlikovanja[6]<prag_odluke) || (mere_razlikovanja[7]<prag_odluke) || (mere_razlikovanja[8]<prag_odluke) )
				{broj_pogresnih_prihvatanja+=1;}*/
				
				//cout<<"Test "<<redni_broj_testa<<" /42732 "<<"  Referentni snimak "<<referentni_snimak<<"  Test snimak "<<test_snimak<<"  Mera razlikovanja "<<mera_razlikovanja<<"  Broj pogresnih prihvatanja "<<broj_pogresnih_prihvatanja<<" / "<<br_testova_pogr_prihv<<" /(41580)"<<endl;}
			cout<<"Test "<<redni_broj_testa<<" /42444 za referentnog govornika "<<referentni_govornik<<"  Test snimak "<<test_snimak<<
				"  Broj pogresnih prihvatanja "<<broj_pogresnih_prihvatanja<<" / "<<br_testova_pogr_prihv<<" /(41580)"<<endl;}
			//cout<<"Test snimak "<<test_snimak<<" Mera razlikovanja "<<mera_razlikovanja<<"   Broj pogresnih prihvatanja "<<broj_pogresnih_prihvatanja<<endl;}
		}//test snimak != referentni snimak
		}//j
	cout << "\n\nTrenutna vrednost verovatnoce pogresnog odbacivanja iznosi " << broj_pogresnih_odbacivanja << " / " << br_testova_pogr_odb
		 << " a verovatnoce pogresnog prihvatanja iznosi " << broj_pogresnih_prihvatanja << " / " << br_testova_pogr_prihv << endl;
	}//i
		p_pogresnog_odbacivanja = (double)broj_pogresnih_odbacivanja/(double)br_testova_pogr_odb;
		p_pogresnog_prihvatanja = (double)broj_pogresnih_prihvatanja/(double)br_testova_pogr_prihv;
		cout<<"\nPrepoznavac: "<<brkropsega<<" e_d kriticnih opsega, s=2, "<<brojMFCCs<<" MFCCs + "<<prvo_dodatno_obelezje<<
			"*e1 + "<<drugo_dodatno_obelezje<<"*e2"<<endl;
		cout << "\nZa vrednost praga odluke " << prag_odluke << "\nVerovatnoca pogresnog odbacivanja iznosi " << p_pogresnog_odbacivanja << " a verovatnoca pogresnog prihvatanja iznosi " << p_pogresnog_prihvatanja << endl;
			//ako pripada i-tom govorniku onda se racuna verovatnoca pogresnog odbacivanja 
			//ako ne pripada i-tom govorniku onda se racuna verovatnoca pogresnog prihvatanja
	
	}//nacin rada = 5  ----  kraj

	if(nacin_rada==6)
	{//u ovom delu bi trebalo obezbediti da program izvrsi testiranje nad govornom bazom izuzimajuci fajlove koji su korisceni pri obuci
		//npr. fajlovi koji u svom nazivu imaju s15, npr. fajl frf01_s15_solo.wav, su korisceni pri obuci za svakog govornika
		//onda bi svi ostali fajlovi trebali biti korisceni pri testiranju
		//string oznaka_referentnog_snimka1 = "s15"; string oznaka_referentnog_snimka2 = "s33"; string oznaka_referentnog_snimka3="s01";
		//string oznaka_referentnog_snimka4 = "s07"; string oznaka_referentnog_snimka5 = "s05"; string oznaka_referentnog_snimka6="s31";
		//string oznaka_referentnog_snimka7 = "s18"; string oznaka_referentnog_snimka8 = "s20"; string oznaka_referentnog_snimka9="s25";
		//string oznaka_referentnog_snimka10 = "s27";string oznaka_referentnog_snimka11 = "s28";string oznaka_referentnog_snimka12 = "s26";
		//string oznaka_referentnog_snimka13 = "s29";string oznaka_referentnog_snimka14 = "s30";string oznaka_referentnog_snimka15 = "s21";
		//string oznaka_referentnog_snimka16 = "s24";

   string oznaka_test_snimka1="29";string oznaka_test_snimka2="30";string oznaka_test_snimka3="31";string oznaka_test_snimka4="32";

		string snimak = "";//ovde ce se upisivati ime snimka koji ce se analizirati
		ifstream snimci;  snimci.open("snimci.txt");//Fajlovi.txt sadrzi nazive svih snimaka u razmatranoj govornoj bazi
		int broj_tacno_prepoznatih = 0;
		int broj_uradjenih_testova = 0;//brojac uradjenih testova, testiranje ce biti uradjeno ako snimak nije koriscen pri obuci 
		for(int i=0;i<11369;i++) //ukupno ima 36*37=1332 snimka u govornoj bazi
		{
			getline(snimci, snimak);//prebacuje se ime snimka iz .txt fajla preko toka (strima) snimci u string snimak
//			string osecanje="";//za ispitivanje kom osecanju pripada snimak
//			string::iterator pocetak_o = snimak.begin()+62;
//			string::iterator kraj_o = snimak.begin()+63;
//			osecanje.append(pocetak_o,kraj_o);
//			if(osecanje=="n") {//L,n,r,S,T
			//sada prvo treba utvrditi da li snimak u svom nazivu ima oznaku referentnog test snimka
	string poredbeni; poredbeni = ""; //radi utvrdjivanja na osnovu imena snimka da li ima ili nema oznaku referentnog snimka
			//te da li se stoga taj snimak posmatra pri testiranju
	string::iterator pocetak = snimak.begin()+60;//6   // http://anaturb.net/C/string_exapm.htm
	string::iterator kraj = snimak.begin()+62;//9
	poredbeni.append(pocetak,kraj);//iz naziva snimka se vadi string izmedju 7. i 9. znaka, to je podstring koji sadrzi 7, 8 i 9. znak
			//cout << poredbeni << endl; //taj podstring se upisuje u string poredbeni
			
			//if ((poredbeni != oznaka_referentnog_snimka1) && (poredbeni != oznaka_referentnog_snimka2) && (poredbeni != oznaka_referentnog_snimka3) && (poredbeni != oznaka_referentnog_snimka4) && (poredbeni != oznaka_referentnog_snimka5) && (poredbeni != oznaka_referentnog_snimka6) && (poredbeni != oznaka_referentnog_snimka7) && (poredbeni != oznaka_referentnog_snimka8) && (poredbeni != oznaka_referentnog_snimka9) && (poredbeni != oznaka_referentnog_snimka10) && (poredbeni != oznaka_referentnog_snimka11) && (poredbeni != oznaka_referentnog_snimka12) && (poredbeni != oznaka_referentnog_snimka13) && (poredbeni != oznaka_referentnog_snimka14) && (poredbeni != oznaka_referentnog_snimka15) && (poredbeni != oznaka_referentnog_snimka16) && (poredbeni != "f01") && (poredbeni != "f02") && (poredbeni != "f03") && (poredbeni != "f04"))
	//if ((poredbeni != oznaka_referentnog_snimka1) && (poredbeni != oznaka_referentnog_snimka2) && (poredbeni != oznaka_referentnog_snimka3) && (poredbeni != oznaka_referentnog_snimka4) && (poredbeni != oznaka_referentnog_snimka5) && (poredbeni != oznaka_referentnog_snimka6) && (poredbeni != oznaka_referentnog_snimka7) && (poredbeni != oznaka_referentnog_snimka8) && (poredbeni != oznaka_referentnog_snimka9) && (poredbeni != "f01") && (poredbeni != "f02") && (poredbeni != "f03") && (poredbeni != "f04"))
				//if ((poredbeni != oznaka_referentnog_snimka1) && (poredbeni != "f01") && (poredbeni != "f02") && (poredbeni != "f03") && (poredbeni != "f04"))
           if((poredbeni==oznaka_test_snimka1)||(poredbeni==oznaka_test_snimka2)||(poredbeni==oznaka_test_snimka3)
			   ||(poredbeni==oznaka_test_snimka4))
			{
				broj_uradjenih_testova++;
				//cout << poredbeni << endl;
				//nad ovim snimcima ce se vrsiti testiranje, jos treba iz imena svakog snimka izdvojiti ime govornika kom pripada
				//i na osnovu poredjenja imena izdvojenog iz naziva snimka i imena govornika koji je prepoznat uvecavati broj tacno prepoznatih
				string govornik = ""; //string za cuvanje imena govornika kom pripada testirani snimak
				string::iterator pocetak1 = snimak.begin()+43;//+0 //pocetak imena govornika u nazivu snimka
				string::iterator kraj1 = snimak.begin()+47;//+5 //zavrsetak imena govornika u nazivu snimka
				govornik.append(pocetak1,kraj1);
				//testiranje:ucitavanje wav snimka, pravljenje vektora obelezja i modela i prepoznavanje kom govorniku pripada
				string prepoznati = ""; //za cuvanje prepoznatog govornika
				if(_access("govornici.txt",0)==-1){cout<<"Baza govornika je prazna, prvo pokrenuti obuku, opcija 1" << endl;}
				else{
					vector<vector<vector<double>>> kovmat1;
					kovmat1.swap(ucitavanje_modela(broj_govornika,broj_obelezja,"modeli.txt"));//u kovmat1 se kopiraju modeli dobijeni pri obuci
					
					vector<double> odbirci_test; vector<double> veksrvr_test; vector<vector<double>> kovmat_test;
		//pravljenje modela test snimka koji se ispituje
					odbirci_test.swap(ucitavanjewavfajla(snimak.c_str()));
					brojprozora = ((odbirci_test.size())-N)/pomerajprozora;
					vector<vector<double>> mfccs_test(brojprozora, vector<double> (brojMFCCs));
					mfccs_test.swap(mfcc(N, pomerajprozora, odbirci_test, oblik_kr_opsega, brkropsega, brojMFCCs, obrnutired, prozor, costab, sintab, costabMFCC));
		vector<vector<double>> vektori_dodatnih_obelezja(brojprozora,vector<double> (3));
		vektori_dodatnih_obelezja.swap(dodatna_obelezja(N, pomerajprozora, odbirci_test, obrnutired, prozor, costab, sintab,prvo_dodatno_obelezje,drugo_dodatno_obelezje,broj_dodatnih_obelezja));
		vector<vector<double>> vektori_obelezja(brojprozora, vector<double> (broj_obelezja));
		if(broj_obelezja == brojMFCCs) {vektori_obelezja.swap(mfccs_test);}
		if(broj_obelezja > brojMFCCs) 
		{
			for (int i=0; i<brojprozora; i++)
			{
				for(int j=0; j<broj_obelezja; j++)
				{
					if(j<brojMFCCs) {vektori_obelezja[i][j] = mfccs_test[i][j];}
					//if(j>=brojMFCCs) {vektori_obelezja[i][j] = vektor_dodatnih_obelezja[i];}
					
					//if(j<brojMFCCs-1) {vektori_obelezja[i][j] = mfccs_test[i][j];}
					//if(j==brojMFCCs-1) {vektori_obelezja[i][j] = vektor_dodatnih_obelezja[i];}
					
					if(j>=brojMFCCs) {
						if(broj_dodatnih_obelezja==1 || broj_dodatnih_obelezja==2 || broj_dodatnih_obelezja==3)
						{vektori_obelezja[i][brojMFCCs]=vektori_dodatnih_obelezja[i][0];}
						if(broj_dodatnih_obelezja==2 || broj_dodatnih_obelezja==3)
						{vektori_obelezja[i][brojMFCCs+1]=vektori_dodatnih_obelezja[i][1];}
						if(broj_dodatnih_obelezja==3)
						{vektori_obelezja[i][brojMFCCs+2]=vektori_dodatnih_obelezja[i][2];}
					}
				}
			}
		}
					veksrvr_test.swap(veksrvr(broj_obelezja,vektori_obelezja));
					kovmat_test.swap(kovmat(broj_obelezja,vektori_obelezja,veksrvr_test));

					vector<double> vekapsodstkovmat_test(broj_govornika);   //vektor razlikovanja modela test govora i modela dobijenih prilikom obuke
	
					for(int i=0; i<broj_govornika;i++)
					{
						vekapsodstkovmat_test[i] = (apsodstupanjekovmat(broj_obelezja, kovmat1[i], kovmat_test));
					}

					prepoznati = prepoznavanje_nad_snimkom_s(snimak, broj_govornika, vektor_govornika, vekapsodstkovmat_test);
					if(prepoznati == govornik) {broj_tacno_prepoznatih++;}
					cout << "Trenutna tacnost prepoznavanja " << broj_tacno_prepoznatih << " / " << broj_uradjenih_testova << endl;
				} //zavrsetak za else
			} //zavrsetak za if(poredbeni==...
//			} //zavrsetak za if(osecanje=="n")
		} //zavrsetak za petlju po i
		cout<<"\nPrepoznavac: "<<brkropsega<<" "<<oblik_kr_opsega<<" kriticna opsega, "<<brojMFCCs<<" MFCCs + "<<prvo_dodatno_obelezje<<
			"*e1 + "<<drugo_dodatno_obelezje<<"*e2 + "<<trece_dodatno_obelezje<<"*e3"<<endl;
		cout << "\nTacnost prepoznavanja " << broj_tacno_prepoznatih << " / " << broj_uradjenih_testova << "  =  " <<
			((double)broj_tacno_prepoznatih/(double)broj_uradjenih_testova)*100.0 << "%" << endl;
	}//nacin_rada=6 --- kraj

	if(nacin_rada==7) //mogu se iskoristiti delovi koda za nacin_rada==3
	{//vrsi se obuka odnosno pravljenje referentnih modela za vektore obelezja sastavljene na osnovu analize energija
		ofstream baza_g;
		if(_access("govornici.txt",0)==-1){baza_g.open("govornici.txt");baza_g<<ime_govornika;baza_g.close();broj_govornika=1;} //ako nema baze govornika i program je prvi put pokrenut
		else{baza_g.open("govornici.txt",std::ios::app);baza_g<<"\n"<<ime_govornika;baza_g.close();broj_govornika++;} //ubacuje ime govornika za koga ce se vrsiti obuka u bazu imena govornika
			
		vector<double> odbirci; vector<double> vek_sr_vr; vector<vector<double>> kov_mat;
		
		odbirci.swap(ucitavanjewavfajla(snimak_za_obuku));
		brojprozora = ((odbirci.size())-N)/pomerajprozora; //racunanje modela snimljenog fajla za obuku
		
		vector<vector<double>> vektori_obelezja(brojprozora,vector<double> (broj_energetskih_obelezja));
		vektori_obelezja.swap(obelezja(N, pomerajprozora, odbirci, obrnutired, prozor, costab, sintab,broj_energetskih_obelezja));
		
		vek_sr_vr.swap(veksrvr(broj_energetskih_obelezja,vektori_obelezja));
		kov_mat.swap(kovmat(broj_energetskih_obelezja,vektori_obelezja,vek_sr_vr));
		
		ofstream modeli; modeli.open("modeli.txt", std::ios::app); // model se dodaje na kraj txt baze modela
		modeli << "\n" << "s" << ime_govornika << "\n";  //ubacivanje imena govornika u modeli.txt fajl iza koga slede elementi kovarijansne matrice, skupa za obuku i test fajla
		for (int i=0; i<broj_energetskih_obelezja; i++)
		{
			for (int j=0; j<broj_energetskih_obelezja; j++)
			{
				modeli << "KovMat " << i+1 << ", " << j+1 << "   " << kov_mat[i][j] << "   " << kov_mat[i][j] << "   " << kov_mat[i][j]-kov_mat[i][j] << "\n";  //radi uporednog prikaza kovmat dobijene pri obuci i za test fajl
			}
		}
		modeli.close();
		cout<<"\nObuka prepoznavaca na govornika '"<<ime_govornika<<"' je uspesno izvrsena i sada u bazi modela ima "<<broj_govornika<< " modela"<<endl;
	}//nacin_rada=7 --- kraj
	if(nacin_rada==8) //mogu se iskoristiti delovi koda za nacin_rada==6
	{//ovde ce se vrsiti provera tacnosti prepoznavanja za vektor obelezja sastavljen na osnovu analize energija u signalu
		//moglo bi se ovde obezbediti da se vrsi i obuka i test posto ce se vrsiti eksperimentisanje nad obelezjima
		//u ovom delu bi trebalo obezbediti da program izvrsi testiranje nad govornom bazom izuzimajuci fajlove koji su korisceni pri obuci
		//npr. fajlovi koji u svom nazivu imaju s15, npr. fajl frf01_s15_solo.wav, su korisceni pri obuci za svakog govornika
		//onda bi svi ostali fajlovi trebali biti korisceni pri testiranju
		string oznaka_referentnog_snimka = "s15";
		string snimak = "";//ovde ce se upisivati ime snimka koji ce se analizirati
		ifstream snimci;  snimci.open("Fajlovi.txt");
		int broj_tacno_prepoznatih = 0;
		int broj_uradjenih_testova = 0;//brojac uradjenih testova, testiranje ce biti uradjeno ako snimak nije koriscen pri obuci 
		for(int i=0;i<1332;i++)
		{
			getline(snimci, snimak);//prebacuje se ime snimka iz .txt fajla preko toka (strima) snimci u string snimak
			//sada prvo treba utvrditi da li snimak u svom nazivu ima oznaku referentnog snimka
			string poredbeni; poredbeni = ""; //radi utvrdjivanja na osnovu imena snimka da li ima ili nema oznaku referentnog snimka
			//te da li se stoga taj snimak posmatra pri testiranju
			string::iterator pocetak = snimak.begin()+6;   // http://anaturb.net/C/string_exapm.htm
			string::iterator kraj = snimak.begin()+9;
			poredbeni.append(pocetak,kraj);//iz naziva snimka se vadi string izmedju 7. i 9. znaka, to je podstring koji sadrzi 7, 8 i 9. znak
			//cout << poredbeni << endl;
			if ((poredbeni != oznaka_referentnog_snimka) && (poredbeni != "f01") && (poredbeni != "f02") && (poredbeni != "f03") &&(poredbeni != "f04"))
			{
				broj_uradjenih_testova++;
				//cout << poredbeni << endl;
				//nad ovim snimcima ce se vrsiti testiranje, jos treba iz imena svakog snimka izdvojiti ime govornika kom pripada
				//i na osnovu poredjenja imena izdvojenog iz naziva snimka i imena govornika koji je prepoznat uvecavati broj tacno prepoznatih
				string govornik = ""; //string za cuvanje imena govornika kom pripada testirani snimak
				string::iterator pocetak1 = snimak.begin(); //pocetak imena govornika u nazivu snimka
				string::iterator kraj1 = snimak.begin()+5; //zavrsetak imena govornika u nazivu snimka
				govornik.append(pocetak1,kraj1);
				//testiranje:ucitavanje wav snimka, pravljenje vektora obelezja i modela i prepoznavanje kom govorniku pripada
				string prepoznati = ""; //za cuvanje prepoznatog govornika
				if(_access("govornici.txt",0)==-1){cout<<"Baza govornika je prazna, prvo pokrenuti obuku, opcija 1" << endl;}
				else{
					vector<vector<vector<double>>> kovmat1;
					kovmat1.swap(ucitavanje_modela(broj_govornika,broj_energetskih_obelezja,"modeli.txt"));
					
					vector<double> odbirci_test; vector<double> veksrvr_test; vector<vector<double>> kovmat_test;
		
					odbirci_test.swap(ucitavanjewavfajla(snimak.c_str()));
					brojprozora = ((odbirci_test.size())-N)/pomerajprozora;
					
		
					vector<vector<double>> vektori_obelezja(brojprozora,vector<double> (broj_energetskih_obelezja));
					vektori_obelezja.swap(obelezja(N, pomerajprozora, odbirci_test, obrnutired, prozor, costab, sintab,broj_energetskih_obelezja));
					
					veksrvr_test.swap(veksrvr(broj_energetskih_obelezja,vektori_obelezja));
					kovmat_test.swap(kovmat(broj_energetskih_obelezja,vektori_obelezja,veksrvr_test));

					vector<double> vekapsodstkovmat_test(broj_govornika);   //vektor razlikovanja modela test govora i modela dobijenih prilikom obuke
	
					for(int i=0; i<broj_govornika;i++)
					{
						vekapsodstkovmat_test[i] = (apsodstupanjekovmat(broj_energetskih_obelezja, kovmat1[i], kovmat_test));
					}

					prepoznati = prepoznavanje_nad_snimkom_s(snimak, broj_govornika, vektor_govornika, vekapsodstkovmat_test);
					if(prepoznati == govornik) {broj_tacno_prepoznatih++;}
					cout << "Trenutna tacnost prepoznavanja " << broj_tacno_prepoznatih << " / " << broj_uradjenih_testova << endl;
				} //zavrsetak za else
			} //zavrsetak za if
		} //zavrsetak za petlju po i
		cout << "Tacnost prepoznavanja " << broj_tacno_prepoznatih << " / " << broj_uradjenih_testova << endl;
	}//nacin_rada=8 --- kraj

	if(nacin_rada==9) //odredjivanje kvadrata modula DFT (energetskog spektra) po prozorima signala koji je u komandnoj liniji unesen kao signal za obuku
	{
		vector<double> odbirci; vector<double> ener_spek;
		
		odbirci.swap(ucitavanjewavfajla(snimak_za_obuku));
		brojprozora = ((odbirci.size())-N)/pomerajprozora; //broj prozora za koje ce se racunati spektar
		vector<double> izdvojeni_odbirci(N); //za odbirke jednog prozora signala

		ofstream e_spektri; e_spektri.open("e_spektri.txt", std::ios::app); // e_spektar se dodaje na kraj txt baze e_spektara
		e_spektri << "\nEnergetski spektri po prozorima za signal " << snimak_za_obuku << "\n";
		
		for(int i=0;i<brojprozora;i++)
		{
			for(int j=0;j<N;j++)
			{
				izdvojeni_odbirci[j]=odbirci[i*pomerajprozora+j];
			}
		ener_spek.swap(kmfft(N, izdvojeni_odbirci, obrnutired, prozor, costab, sintab));

		//deo koji odredjuje maksimum u celom spektru prozora
		double maksimum = 0.0; //maksimum u prozoru se odredjuje da bi se vrednosti enrgetskog spektra mogle iskazati i u decibelima [dB], pri cemu je maksimalna vrednost u prozoru uzeta za referentnu pri racunanju energetskog spektra u dB
		for(int j=0;j<N/2;j++)  { if(ener_spek[j]>maksimum) {maksimum=ener_spek[j];} }

		//deo koji odredjuje maksimum1 i njegovu ucestanost k1
		int k1=0; double maksimum1 = -10000.0;
		for(int j=25;j<N/2;j++)
		{
			if(log(ener_spek[j])>maksimum1)
			{
				maksimum1 = log(ener_spek[j]); k1 = j;
			}
		}

		//deo koji odredjuje maksimum2 i njegovu ucestanost k2
		int k2 = 0;
		double maksimum2 = -10000;
		for(int j=25;j<N/2;j++) //racunanje drugog po redu maksimuma u spektru
		{
			if(j>k1+10)
			{
				if(log(ener_spek[j])>maksimum2)
				{
					maksimum2 = log(ener_spek[j]); k2=j;
				}
			}
		}

		e_spektri << "\n" << i << ". prozor:" << "\n";  //ubacivanje rednog broja prozora u e_spektri.txt fajl iza koga slede elementi kvadrata modula DFT
		
		for (int j=0; j<N; j++)
		{
			double ener_spek_db = 10.0*log10(ener_spek[j]/maksimum);
			if(j==k1) {e_spektri << "k=k1= " << j << ":   " << ener_spek[j] << "   " << ener_spek_db << " dB" << "\n";}
			if(j==k2) {e_spektri << "k=k2= " << j << ":   " << ener_spek[j] << "   " << ener_spek_db << " dB" << "\n";}
			e_spektri << "k= " << j << ":   " << ener_spek[j] << "   " << ener_spek_db << " dB" << "\n";
		}
		} //kraj petlje po i
		
		e_spektri.close();
		cout<<"\nUpis enrgetskih spektara signala je zavrsen"<<endl;
	}//kraj za if(nacin_rada==9)

	if(nacin_rada==10) //ovde ce se racunati promenljivost obelezja za govornika i izmedju govornika
		//za vektor obelezja koji sadrzi MFCCs+e1+e2
	{   //izracunava se srednje ukupno apsolutno razlikovanje na nivou obelezja za snimak koji se unosi kao snimak za obuku, snimak_za_obuku
		//i za snimak koji se unosi kao test snimak, snimak_za_prepoznavanje
		
		//racunanje obelezja za prvi snimak
		vector<double> odbirci1;
		
		odbirci1.swap(ucitavanjewavfajla(snimak_za_obuku));//ucitavanje prvog snimka, unetog kao snimak_za_obuku
		brojprozora = ((odbirci1.size())-N)/pomerajprozora; 
		
		vector<vector<double>> mfccs1(brojprozora, vector<double> (brojMFCCs));//racunanje vektora MFCCs
		mfccs1.swap(mfcc(N, pomerajprozora, odbirci1, oblik_kr_opsega, brkropsega, brojMFCCs, obrnutired, prozor, costab, sintab, costabMFCC));
		
		vector<vector<double>> vektori_dodatnih_obelezja1(brojprozora,vector<double> (broj_dodatnih_obelezja));//racunanje vektora dodatnih obelezja
		vektori_dodatnih_obelezja1.swap(dodatna_obelezja(N, pomerajprozora, odbirci1, obrnutired, prozor, costab, sintab,prvo_dodatno_obelezje,drugo_dodatno_obelezje,broj_dodatnih_obelezja));	
		
		vector<vector<double>> vektori_obelezja1(brojprozora, vector<double> (broj_obelezja));//spajanje vektora MFCCs i vektora dodatnih obelezja u jedan vektor obelezja
		if(broj_dodatnih_obelezja == 0) {vektori_obelezja1.swap(mfccs1);}
		if(broj_dodatnih_obelezja != 0) 
		{
			for (int i=0; i<brojprozora; i++)
			{
				for(int j=0; j<broj_obelezja; j++)
				{
					if(j<brojMFCCs) {vektori_obelezja1[i][j] = mfccs1[i][j];}
					
					if(j>=brojMFCCs) {
						if(broj_dodatnih_obelezja==1 || broj_dodatnih_obelezja==2 || broj_dodatnih_obelezja==3)
						{vektori_obelezja1[i][brojMFCCs]=vektori_dodatnih_obelezja1[i][0];}
						if(broj_dodatnih_obelezja==2 || broj_dodatnih_obelezja==3)
						{vektori_obelezja1[i][brojMFCCs+1]=vektori_dodatnih_obelezja1[i][1];}
						if(broj_dodatnih_obelezja==3)
						{vektori_obelezja1[i][brojMFCCs+2]=vektori_dodatnih_obelezja1[i][2];}
					}
				}
			}
		}
		//zavrsetak racunanja obelezja za prvi snimak unetog kao snimak_za_obuku

		int brojprozora1=brojprozora;//radi cuvanja broja prozora za prvi snimak, pre nego sto pocne racunanje obelezja za drugi snimak
		//u prethodnu vrednost promenljive brojprozora bice upisana nova vrednost za drugi snimak
		//ovo je ostavljeno ovako posto se promenljiva broj prozora kao takva koristi pri racunanju obelezja
		//posto je brojprozora na pocetku ovog celog programa definisana kao globalna promenljiva

		//racunanje obelezja za drugi snimak
		vector<double> odbirci2;
		
		odbirci2.swap(ucitavanjewavfajla(snimak_za_prepoznavanje));//ucitavanje drugog snimka, unetog kao snimak_za_prepoznavanje
		brojprozora = ((odbirci2.size())-N)/pomerajprozora; 
		
		vector<vector<double>> mfccs2(brojprozora, vector<double> (brojMFCCs));//racunanje vektora MFCCs
		mfccs2.swap(mfcc(N, pomerajprozora, odbirci2, oblik_kr_opsega, brkropsega, brojMFCCs, obrnutired, prozor, costab, sintab, costabMFCC));
		
		vector<vector<double>> vektori_dodatnih_obelezja2(brojprozora,vector<double> (broj_dodatnih_obelezja));//racunanje vektora dodatnih obelezja
		vektori_dodatnih_obelezja2.swap(dodatna_obelezja(N, pomerajprozora, odbirci2, obrnutired, prozor, costab, sintab,prvo_dodatno_obelezje,drugo_dodatno_obelezje,broj_dodatnih_obelezja));	
		
		vector<vector<double>> vektori_obelezja2(brojprozora, vector<double> (broj_obelezja));//spajanje vektora MFCCs i vektora dodatnih obelezja u jedan vektor obelezja
		if(broj_dodatnih_obelezja == 0) {vektori_obelezja2.swap(mfccs2);}
		if(broj_dodatnih_obelezja != 0) 
		{
			for (int i=0; i<brojprozora; i++)
			{
				for(int j=0; j<broj_obelezja; j++)
				{
					if(j<brojMFCCs) {vektori_obelezja2[i][j] = mfccs2[i][j];}
					
					if(j>=brojMFCCs) {
						if(broj_dodatnih_obelezja==1 || broj_dodatnih_obelezja==2 || broj_dodatnih_obelezja==3)
						{vektori_obelezja2[i][brojMFCCs]=vektori_dodatnih_obelezja2[i][0];}
						if(broj_dodatnih_obelezja==2 || broj_dodatnih_obelezja==3)
						{vektori_obelezja2[i][brojMFCCs+1]=vektori_dodatnih_obelezja2[i][1];}
						if(broj_dodatnih_obelezja==3)
						{vektori_obelezja2[i][brojMFCCs+2]=vektori_dodatnih_obelezja2[i][2];}
					}
				}
			}
		}
		//zavrsetak racunanja obelezja za drugi snimak unetog kao snimak_za_prepoznavanje

		int brojprozora2=brojprozora;

		//da bi se izvrsilo poredjenje vektora obelezja za prvi i drugi snimak prethodno cu izjednaciti broj vektora obelezja za prvi i drugi snimak
		//porede se brojprozora1 i brojprozora2 i kao brojprozora ce se uzeti manji od ta dva broja, odnosno isti ako su isti
		if(brojprozora1<=brojprozora2) {brojprozora=brojprozora1;}
		if(brojprozora2<=brojprozora1) {brojprozora=brojprozora2;}

		//sledi deo u kome se racuna razlikovanje
		vector<double> razlikovanje(broj_obelezja); for(int i=0;i<broj_obelezja;i++){razlikovanje[i]=0.0;}

		for(int i=0;i<brojprozora;i++)
		{
			for(int j=0;j<broj_obelezja;j++)
			{
				razlikovanje[j]+=abs(vektori_obelezja2[i][j]-vektori_obelezja1[i][j]);
			}
		}

		for(int i=0;i<broj_obelezja;i++) { //normalizovanje razlikovanja brojem prozora, usrednjavanje u sustini
			razlikovanje[i]=razlikovanje[i]/((double)brojprozora);}
		
		cout<<"\nUsrednjeni zbir apsolutnih razlikovanja po obelezjima:"<<endl;
		for(int i=0;i<broj_obelezja;i++)
		{
			cout<<i+1<<".  "<<razlikovanje[i]<<endl;
		}
	}//zavrsetak za nacin_rada=10

		if(nacin_rada==11)
	{//u ovom delu bi trebalo obezbediti da program izvrsi obuku nad govornom bazom izuzimajuci fajlove koji ce se koristiti pri testu
		//snimci koji u svom nazivu imaju oznake 16n, 17n i 18n su za testiranje
		//onda bi svi ostali snimci trebali biti korisceni pri obuci
		ofstream baza_g;//izlazni tok (strim) podataka za upis u bazu govornika
		ofstream modeli;//izlazni tok (strim) podataka za upis u bazu modela
		vector<double> odbirci; vector<double> vek_sr_vr; vector<vector<double>> kov_mat;
 string oznaka_test_snimka1 = "29";string oznaka_test_snimka2 = "30";string oznaka_test_snimka3="31";string oznaka_test_snimka4="32";

		string snimak = "";//ovde ce se upisivati ime snimka koji ce se analizirati
		ifstream snimci;  snimci.open("snimci.txt");//snimci.txt sadrzi nazive svih snimaka u razmatranoj govornoj bazi
		//ulazni tok (strim) podataka od snimka ka programu
		int broj_modela = 0;//??
		int broj_govornika = 0;
		for(int i=0;i<11369;i++) //ukupno ima 36*37=1332 snimka u govornoj bazi
		{
			getline(snimci, snimak);//prebacuje se ime snimka iz .txt fajla preko toka (strima) snimci u string snimak
			//sada prvo treba utvrditi da li snimak u svom nazivu ima oznaku test snimka
			string poredbeni; poredbeni = ""; //radi utvrdjivanja na osnovu imena snimka da li ima ili nema oznaku test snimka
			//te da li se stoga taj snimak posmatra pri obuci
			string::iterator pocetak = snimak.begin()+60;   // http://anaturb.net/C/string_exapm.htm
			string::iterator kraj = snimak.begin()+62;
			poredbeni.append(pocetak,kraj);//iz naziva snimka se vadi string izmedju 61. i 62. znaka,
			               //to je podstring koji sadrzi 61. i 62. znak
			//cout << poredbeni << endl; //taj podstring se upisuje u string poredbeni
//			string oznaka; oznaka="";//za vrstu osecanja(n,L,r,S,T)
//			string::iterator pocetak1=snimak.begin()+62; string::iterator kraj1=snimak.begin()+63;
//			oznaka.append(pocetak1,kraj1);
			
			if ((poredbeni != oznaka_test_snimka1) && (poredbeni != oznaka_test_snimka2) && (poredbeni != oznaka_test_snimka3)
				&& (poredbeni != oznaka_test_snimka4))// && (oznaka == "n"))//tada se vrsi obuka
			{
				broj_modela++;
				//cout << poredbeni << endl;
				//ovi snimci su za obuku, jos treba iz imena svakog snimka izdvojiti ime govornika kom pripada
				string govornik = ""; //string za cuvanje imena govornika kom pripada snimak za obuku
				string::iterator pocetak2 = snimak.begin()+43; //pocetak imena govornika u nazivu snimka
				string::iterator kraj2 = snimak.begin()+47; //zavrsetak imena govornika u nazivu snimka
				govornik.append(pocetak2,kraj2);
				//ofstream baza_g; izlazni tok za upis u bazu govornika govornici.txt
		if(_access("govornici.txt",0)==-1){baza_g.open("govornici.txt");baza_g<<govornik;baza_g.close();broj_govornika=1;} //ako nema baze govornika i program je prvi put pokrenut
		else{baza_g.open("govornici.txt",std::ios::app);baza_g<<"\n"<<govornik;baza_g.close();broj_govornika++;} //ubacuje ime govornika za koga ce se vrsiti obuka u bazu imena govornika
				cout<<"\nPravljenje modela za govornika "<<govornik<<" na osnovu snimka "<<snimak<<endl;
				//obuka:ucitavanje wav snimka, pravljenje vektora obelezja i modela i upisivanje u bazu modela govornika
					
				//vector<double> odbirci; vector<double> veksrvr; vector<vector<double>> kovmat;
		//pravljenje modela snimka koji se ispituje
				odbirci.swap(ucitavanjewavfajla(snimak.c_str()));
				brojprozora = ((odbirci.size())-N)/pomerajprozora;
				vector<vector<double>> mfccs(brojprozora, vector<double> (brojMFCCs));
				mfccs.swap(mfcc(N, pomerajprozora, odbirci, oblik_kr_opsega, brkropsega, brojMFCCs, obrnutired, prozor, costab, sintab, costabMFCC));
		        vector<vector<double>> vektori_dodatnih_obelezja(brojprozora,vector<double> (3));
		        vektori_dodatnih_obelezja.swap(dodatna_obelezja(N, pomerajprozora, odbirci, obrnutired, prozor, costab, sintab,prvo_dodatno_obelezje,drugo_dodatno_obelezje,broj_dodatnih_obelezja));
		        vector<vector<double>> vektori_obelezja(brojprozora, vector<double> (broj_obelezja));
		        if(broj_obelezja == brojMFCCs) {vektori_obelezja.swap(mfccs);}
		if(broj_obelezja > brojMFCCs) 
		{
			for (int i=0; i<brojprozora; i++)
			{
				for(int j=0; j<broj_obelezja; j++)
				{
					if(j<brojMFCCs) {vektori_obelezja[i][j] = mfccs[i][j];}
					
					if(j>=brojMFCCs) {
						if(broj_dodatnih_obelezja==1 || broj_dodatnih_obelezja==2 || broj_dodatnih_obelezja==3)
						{vektori_obelezja[i][brojMFCCs]=vektori_dodatnih_obelezja[i][0];}
						if(broj_dodatnih_obelezja==2 || broj_dodatnih_obelezja==3)
						{vektori_obelezja[i][brojMFCCs+1]=vektori_dodatnih_obelezja[i][1];}
						if(broj_dodatnih_obelezja==3)
						{vektori_obelezja[i][brojMFCCs+2]=vektori_dodatnih_obelezja[i][2];}
					}
				}
			}
		}
					vek_sr_vr.swap(veksrvr(broj_obelezja,vektori_obelezja));
					kov_mat.swap(kovmat(broj_obelezja,vektori_obelezja,vek_sr_vr));
					//upisivanje imena govornika i modela u bazu modela govornika (u modeli.txt)
					modeli.open("modeli.txt", std::ios::app); // model se dodaje na kraj txt baze modela
		modeli << "\n" << "s" << govornik << "\n";  //ubacivanje imena govornika u modeli.txt fajl iza koga slede elementi kovarijansne matrice, skupa za obuku i test fajla
		for (int i=0; i<broj_obelezja; i++)
		{
			for (int j=0; j<broj_obelezja; j++)
			{
				modeli << "KovMat " << i+1 << ", " << j+1 << "   " << kov_mat[i][j] << "   " << kov_mat[i][j] << "   " << kov_mat[i][j]-kov_mat[i][j] << "\n";  //radi uporednog prikaza kovmat dobijene pri obuci i za test fajl
			}
		}
		modeli.close();
		cout<<"\nObuka prepoznavaca na govornika '"<<govornik<<"' je uspesno izvrsena i sada u bazi modela ima "<<broj_govornika<< " modela"<<endl;
			} //zavrsetak za if
		} //zavrsetak za petlju po i
	}//nacin_rada=11 --- kraj

	//cout << "\nDa bi zatvorili prozor programa pritisnite taster Enter. \n";
	//cin.ignore ( numeric_limits<streamsize>::max(), '\n' ); //Clean the stream and ask for input
	//cin.get();  //kad bi se samo cin.get() koristilo bio bi potreban samo jedan Enter za zatvaranje prozora aplikacije
                //ukoliko je vec ranije u programu vrsen unos preko tastature, moze biti da je pri unosu koriscen i Enter te on onda vec postoji u strimu
                //zato se primenjuje brisanje ulaznog strima pomocu cin.ignore
                //prozor programa se zatvara nakon dva pritiska tastera Enter
                //preuzeto sa:  http://www.dreamincode.net/forums/topic/30581-holding-the-execution-window-open/



	return 0;
} 


/* #include <mmsystem.h>
UINT waveInGetNumDevs(VOID); 
http://msdn.microsoft.com/en-us/library/ms177549(v=VS.80).aspx */

//http://stackoverflow.com/questions/644673/fast-way-to-copy-one-vector-into-another*/      //.swap

//unos eha...
		/*vector<float> zakasnjenje(4410);   //2205 za 100ms
		for (int j=0; j<4410; j++)
		{
			zakasnjenje[j] = 0.0;
		}
		vector<float> zakasnjenje1;
		zakasnjenje1.swap(zakasnjenje);
		zakasnjenje1.insert(zakasnjenje1.end(), odbirci1.begin(), odbirci1.end());
		odbirci1.insert(odbirci1.end(), zakasnjenje.begin(), zakasnjenje.end());
		int broj;
		broj = odbirci1.size();
		for (int j=0; j<broj; j++)
		{
			odbirci1[j] = odbirci1[j] + zakasnjenje1[j];
		}*/
		//unos eha.