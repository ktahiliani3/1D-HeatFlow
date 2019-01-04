#include<iostream>
#include<fstream>
#include<iomanip>
#include<limits>
#include"mpi.h"
#include<math.h>

using namespace std;


int main (int argc, char**argv){

    int rc, processors, rank;
    rc = MPI_Init(&argc, &argv);
    double T1, T2, r, k = 1, h = 2;


    T1 = atof(argv[1]);
    T2 = atof(argv[2]);
    int GridPoints = atoi(argv[3]);
    int TimeSteps = atoi(argv[4]);



    r = k/(h*h);

    if (rc !=MPI_SUCCESS){

        cout<<"Unable to start MPI";
        MPI_Abort(MPI_COMM_WORLD, rc);

    }
    else{

        int ProcessorStart[processors];
        int ProcessorEnd[processors];
        int NumberofPoints[processors];



        MPI_Comm_size(MPI_COMM_WORLD, &processors);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        int SlavePoints = GridPoints/processors;
        int MasterPoints = GridPoints/processors + GridPoints%processors;


        if (processors<=GridPoints){



            for (int i = 0; i < processors - 1; i++){

                ProcessorStart[i] = i*SlavePoints;
                ProcessorEnd[i] = ProcessorStart[i] + SlavePoints - 1;
                NumberofPoints[i] = SlavePoints + 2;

            }

            ProcessorStart[processors-1] = ProcessorStart[processors-2] + SlavePoints;
            ProcessorEnd[processors-1] = GridPoints - 1;
            NumberofPoints[processors-1] = MasterPoints + 2;
        }
        else{

           for(int i = 0; i < processors; i++){
                if(i<GridPoints){
                    ProcessorStart[i] = i;
                    ProcessorEnd[i] = i+1;
                    NumberofPoints[i] = SlavePoints + 3;
                }
                else{
                    ProcessorStart[i] = i;
                    ProcessorEnd[i] = i+1;
                    NumberofPoints[i] = 0;
                }
           }

           processors = GridPoints;

        }
        GridPoints = GridPoints + 2;

        int CurrentWidth = NumberofPoints[rank];
        int CurrentStart = ProcessorStart[rank];
        int CurrentEnd = ProcessorEnd[rank];

        double CurrentGrid[CurrentWidth];
        double CurrentTemp[CurrentWidth];

        if (rank == 0)
        {
            CurrentGrid[0] = T1;

            if (processors==1){
                CurrentGrid[CurrentWidth -1] = T2;
            }
        }
        else{


            if(rank == processors - 1){
                CurrentGrid[CurrentWidth - 1] = T2;
            }


        }

        for(int i = 0;i < TimeSteps; i++){
            for(int j = 1; j< CurrentWidth-1;j++){

                CurrentTemp[j] = (1-2*r)*CurrentGrid[j] + r*CurrentGrid[j-1] + r*CurrentGrid[j+1];

            }

            if(((rank % 2) == 0) && (rank < processors)){

                if (rank < processors -1){
                    MPI_Send(&CurrentTemp[CurrentWidth - 2],1, MPI_DOUBLE,rank+1,2,MPI_COMM_WORLD);

                }
                if (rank > 0){
                    MPI_Send(&CurrentTemp[1],1, MPI_DOUBLE,rank-1,2,MPI_COMM_WORLD);

                }
                if (rank < processors - 1){

                    MPI_Status status;

                    MPI_Recv(&CurrentTemp[CurrentWidth-1],1, MPI_DOUBLE,rank+1,2,MPI_COMM_WORLD, &status);

                }
                if(rank > 0){

                    MPI_Status status;

                    MPI_Recv(&CurrentTemp[0],1, MPI_DOUBLE,rank-1,2,MPI_COMM_WORLD, &status);
                }

            }
            else{

                if(rank < processors){
                    if (rank < processors-1){

                        MPI_Status status;

                        MPI_Recv(&CurrentTemp[CurrentWidth-1],1, MPI_DOUBLE,rank+1,2,MPI_COMM_WORLD, &status);

                    }
                    if(rank > 0){

                        MPI_Status status;

                        MPI_Recv(&CurrentTemp[0],1, MPI_DOUBLE,rank-1,2,MPI_COMM_WORLD, &status);
                    }
                    if (rank < processors -1){

                        MPI_Send(&CurrentTemp[CurrentWidth - 2],1, MPI_DOUBLE,rank+1,2,MPI_COMM_WORLD);

                    }
                    if (rank > 0){

                        MPI_Send(&CurrentTemp[1],1, MPI_DOUBLE,rank-1,2,MPI_COMM_WORLD);

                    }
                }
            }
            if (rank == 0){

                CurrentTemp[0] = T1;
            }
            if (rank == processors - 1){

                CurrentTemp[CurrentWidth - 1] = T2;

            }

            MPI_Barrier(MPI_COMM_WORLD);

            for (int p = 0; p < CurrentWidth; p++){

                CurrentGrid[p] = CurrentTemp[p];
            }
        }

        if (rank ==0){

            double NewGrid[GridPoints];

            NewGrid[0] = T1;
            NewGrid[GridPoints - 1] = T2;

            for(int i = 1; i <= CurrentWidth - 2; i++){

                NewGrid[i] = CurrentGrid[i];
            }

            for(int i = 1; i <= processors - 1; i++){

                MPI_Status status;

                MPI_Recv(&NewGrid[0] + ProcessorStart[i]+1,NumberofPoints[i] - 2,MPI_DOUBLE,i, 3,MPI_COMM_WORLD, &status);
            }

            ofstream values ("heat1Doutput.csv", std::ofstream::out);
            for (int i = 0; i < GridPoints - 1; i++){

                values<<NewGrid[i]<<", ";

            }
            values<<NewGrid[GridPoints - 1];
            values.close();


        }
        else{

            if(rank<processors){
                int points = CurrentWidth -2;
                MPI_Send(&CurrentGrid[1],points,MPI_DOUBLE,0,3,MPI_COMM_WORLD);
            }
        }
    }



    MPI_Finalize();
    return EXIT_SUCCESS;
}
























