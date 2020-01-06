
#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>
#include <vector>
#include <ctime>
#include <mpi.h>
#include <math.h>


int const globalBufferLength = 50;

struct distrOpt {
    int gszI, gszJ;         /* Global Size */
    int pnI, pnJ;          /* dim of 2d proc mesh */
    int nIter;             /* Number of iterations */
    int file_jump;         /* Number of iterations before new file gen*/
    int seed;
    int save_file;        /* save to file?*/

};
struct neighbours {
    int left;
    int right;
    int down;
    int up;
};
struct borders {
    std::vector<bool> left;
    std::vector<bool> right;
    std::vector<bool> down;
    std::vector<bool> up;
};


void distr_borders(std::vector<std::vector<bool>> board, neighbours nbr, MPI_Comm Communicator, borders borders);

void initializeBoard(std::vector<std::vector<bool>> &board) {
    int deadCellMultiplyer = 2;
    srand(time(0));
    for (auto &col : board) {
        for (auto element : col) {
            element = (rand() % (deadCellMultiplyer + 1) == 0);
        }
    }
}

void updateBoard(std::vector<std::vector<bool>> &board) {
    const size_t rows = board.size();
    const size_t cols = board[0].size();
    std::vector<std::vector<int>> liveNeighbors(rows, std::vector<int>(cols, 0));

    //Count live neigbors
    for (size_t i=0; i<rows; ++i) {
        for (size_t j=0; j<cols; ++j) {
            if (board[i][j]) {
                for (int di=-1; di<=1; ++di) {
                    for (int dj=-1; dj<=1; ++dj) {
                        //Periodic boundary conditions
                        liveNeighbors[(i+di+rows)%rows][(j+dj+cols)%cols]++;
                    }
                }
                liveNeighbors[i][j]--; //Correction so that a cell does not concider itself as a live neighbor
            }
        }
    }

    //Update board
    for (size_t i=0; i<rows; ++i) {
        for (size_t j=0; j<cols; ++j) {
            board[i][j] = ( (liveNeighbors[i][j] == 3) || (board[i][j] && liveNeighbors[i][j] == 2) );
        }
    }
}

void writeBoardToFile(std::vector<std::vector<bool>> &board, size_t firstRow, size_t lastRow, size_t firstCol, size_t lastCol, std::string fileName, int iteration, uint processID) {
    //Open file
    std::ofstream outputFile(fileName + "_" + std::to_string(iteration) + "_" + std::to_string(processID) + ".gol");
    //Write metadata
    outputFile << std::to_string(firstRow) << " " << std::to_string(lastRow) << std::endl;
    outputFile << std::to_string(firstCol) << " " << std::to_string(lastCol) << std::endl;
    //Write data
    std::ostream_iterator<bool> outputIterator(outputFile, "\t");
    for (size_t i=0; i<board.size(); ++i) {
        copy(board[i].begin(), board[i].end(), outputIterator);
        outputFile << std::endl;
    }
    //Close file
    outputFile.close();
}

std::string setUpProgram(size_t rows, size_t cols, int iteration_gap, int iterations, int processes) {
    //Generate progam name based on current time, all threads should use the same name!
    time_t rawtime;
    struct tm *timeInfo;
    char buffer[globalBufferLength];
    time (&rawtime);
    timeInfo = localtime(&rawtime);
    strftime(buffer, sizeof(buffer), "%Y-%m-%d-%H-%M-%S", timeInfo);
    std::string programName(buffer);

    //Generate main file
    std::ofstream outputFile(programName + ".gol");
    outputFile << std::to_string(rows) << " " << std::to_string(cols) << " " << std::to_string(iteration_gap) << " " << std::to_string(iterations) << " " << std::to_string(processes) << std::endl;
    outputFile.close();

    return programName;
}

int main(int argc, char* argv[]) {
    int  rank, processors;
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &processors);

    MPI_Datatype GOL_options;
    MPI_Aint displ[1] = {sizeof(int)};
    int lengths[1] = {3};
    MPI_Datatype types[1]={MPI_INT};

    MPI_Type_create_struct(3, lengths, displ, types, &GOL_options);
    MPI_Type_commit(&GOL_options);

    distrOpt options;



    if (rank == 0) {
        if (argc != 5) {
            std::cout
                    << "This program should be called with four arguments! \nThese should be, the total number of rows; the total number of columns; the gap between saved iterations and the total number of iterations, in that order."
                    << std::endl;
            MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
            return 1;
        }
        int rows, cols;
        int iteration_gap, iterations;
        try {
            rows = atoi(argv[1]);
            cols = atoi(argv[2]);
            iteration_gap = atoi(argv[3]);
            iterations = atoi(argv[4]);
        } catch (std::exception const &exc) {
            std::cout << "One or more program arguments are invalid!" << std::endl;
            MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
            return 1;
        }
        float z = sqrt(processors);
        if ( rows <= 0 || rows != cols || z !=floor(z )|| rows % (int) z != 0 || rows / (int) z < 4 ){
            std::cout << "Illegal board size parameter combination!" << std::endl;
            MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
            return 1;
        }
        int mesh_dim = (int) z;

        int seed = rand();
         options = {
                cols, rows,            /* Global Size */
                mesh_dim, mesh_dim,          /* dim of 2d proc mesh */
                iterations,             /* Number of iterations */
                iteration_gap,         /* Number of iterations before new file gen*/
                seed,
                0       /* save to file?*/

        };
//        int processes = 1, processID = 0;
//        size_t firstRow = 0, lastRow = rows - 1, firstCol = 0, lastCol = cols - 1;
//        std::string programName = setUpProgram(rows, cols, iteration_gap, iterations, processes);


    }
    MPI_Bcast(&options,1,GOL_options,0,MPI_COMM_WORLD);
    MPI_Type_free(&GOL_options);
    int pcoory = rank / (options.pnI-1);
    int pcoorx = rank % (options.pnJ-1);


    int p_up = pcoory +1;
    int p_down = pcoory -1;
    int p_right = pcoorx +1;
    int p_left = pcoory - 1;

    if (pcoorx == options.pnI-1){
        p_up = -1;
    }
    if (pcoory == options.pnJ-2){
        p_right = -1;
    }

    int left_c = pcoorx*((options.gszI-1)/options.pnI);
    int right_c = (pcoorx + 1)*((options.gszI-1)/options.pnI);
    int down_c = pcoory*((options.gszJ-1)/options.pnJ);
    int up_c = (pcoory+1)*((options.gszJ-1)/options.pnJ);

    neighbours nbr = { p_left,p_right,p_down,p_up};

    //Build board
    std::vector<std::vector<bool>> board(((options.gszI-1)/options.pnI)+3, std::vector<bool>(((options.gszJ-1)/options.pnJ)+3));
    initializeBoard(board);
    borders borders;
    distr_borders(board,nbr,MPI_COMM_WORLD,borders);

    //Do iteration
    writeBoardToFile(board,  firstRow, lastRow, firstCol, lastCol,  programName, 0, processID);
    for (int i=1; i<=iterations; ++i) {
        updateBoard(board);
        if (i%iteration_gap == 0) {
            writeBoardToFile(board,  firstRow, lastRow, firstCol, lastCol,  programName, i, processID);
        }
    }



}

void distr_borders(std::vector<std::vector<bool>> &board, neighbours nbr, MPI_Comm communicator,borders &borders,distrOpt &options){



    static MPI_Datatype mpi_vect = MPI_DATATYPE_NULL;

    int lIsz =(options.gszI-1)/options.pnI;
    int lJsz =(options.gszJ-1)/options.pnJ;

    if (mpi_vect == MPI_DATATYPE_NULL) {
        MPI_Type_vector(lIsz, 1, lJsz+2, MPI_CXX_BOOL, &mpi_vect);
        MPI_Type_commit(&mpi_vect);
    }
    MPI_Request req[4];
    MPI_Isend(board[1][1], 1, mpi_vect ,nbr.left, 0, communicator, req);
    MPI_Isend(board[1][lJsz], 1, mpi_vect ,nbr.right, 0, communicator, req+1);

    MPI_Irecv(board[1][0],1,mpi_vect,nbr.left,0,communicator,req+2);
    MPI_Irecv(board[1][lIsz+1],1,mpi_vect,nbr.left,0,communicator,req+2);

    MPI_Waitall(4,req,MPI_STATUSES_IGNORE);

    MPI_Isend(board[1][0],)



}
