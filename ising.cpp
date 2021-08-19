# include <cmath>
# include <cstdlib>
# include <fstream>
# include <random>
# include <iostream> 



using namespace std;


void create_lattice(int N);
void ising_step(int N, int iterations, float T);
int red[100][100];




int main(){

    int N=30;



    create_lattice(N);
    ising_step(N,1E7, 0.5);
    ising_step(N,1E7, 1);
    ising_step(N,1E7, 1.5);
    ising_step(N,1E7, 2);
    ising_step(N,1E7, 3);
    ising_step(N,1E7, 4);
    ising_step(N,1E7, 6);
    

    N = 100;

    create_lattice(N);
    ising_step(N,1E8, 0.5);
    ising_step(N,1E8, 1);
    ising_step(N,1E7, 1.5);
    ising_step(N,1E7, 2);
    ising_step(N,1E7, 3);
    ising_step(N,1E8, 4);
    ising_step(N,1E8, 6);



}

void create_lattice(int N){
       
   
    random_device rd;
    mt19937 e{rd()};

// Genero la red de spines aleatoriamente, eligiendo entre +-1


    for (int i=0; i<N;i++){
        for (int j=0; j<N;j++){
            uniform_int_distribution<int> dist{0,1};
            red[i][j] = (dist(e)==0 ? -1 : 1);
        }
    }

}

void ising_step(int N, int iterations, float T){

    int ranx, rany, spin, nbrSum;
    double deltaE, beta = 1/T, energia = 0, prob; 
    int mag, J=1, h=1;
    
    random_device rd;
    mt19937 e{rd()};

    ofstream Efile; 
    Efile.open("energia_T="+ to_string(T) +"_N=" + to_string(N) + ".txt"); 


    // Propongo un spin para el flip, calculo el delta de energía y elijo si cambiarlo o no

    ofstream Mfile; 
    Mfile.open("mag_T="+ to_string(T) +"_N=" + to_string(N) + ".txt");
    mag=0;

    for (int i = 0; i<N; i++){
        for (int j = 0; j<N; j++){
            mag += red[i][j];
        }
    }

    for (int i = 0; i<N; i++){
        for (int j = 0; j<N; j++){
            energia  -= J * (red[i][j]*(red[(i+1)%N][j] + red[(i-1+N)%N][j] + red[i][(j+1)%N] + red[i][(j-1+N)%N]));
            energia -= h * mag;
        }
    }
    Efile << "Energías a T=" << T << ", beta="<<beta <<"\n";
    Efile << energia << endl; 

    Mfile << "Magnetizacion a T=" << T << ", N=" << N << "\n";
    Mfile << mag << endl; 

    int deltaMag;

    for(int i = 0 ; i<iterations; i++){

        uniform_int_distribution<int> dist{0,N-1};
        ranx = dist(e);
        rany = dist(e);

        spin=red[ranx][rany];
        nbrSum =  red[(ranx+1)%N][rany] + red[(ranx-1+N)%N][rany] + red[ranx][(rany+1)%N] + red[ranx][(rany-1+N)%N];
        deltaE = 2*spin*nbrSum;
        deltaMag = ((spin>0) - (spin<0))*2;


        if (deltaE<=0){
            red[ranx][rany] *= -1;
            mag -= deltaMag;
            energia += deltaE - h*deltaMag;
        } else
        {
            uniform_real_distribution<float> distr{0,1};
            prob=distr(e);

            if(exp(-deltaE*beta)>prob){
                red[ranx][rany] *= -1;
                mag -= deltaMag;
                energia += deltaE - h*deltaMag;
            }
        }
        Efile << energia << "\n";
        Mfile << mag << "\n";

    }

    cout << "T" << T << " terminado" << "\n";

}
