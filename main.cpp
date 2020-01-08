
#include <iostream>
#include <fstream>
#include <iterator>
#include <vector>
#include <ctime>
#include <mpi.h>
#include <cmath>
#include <string>
#include <thread>

int const globalBufferLength = 50;

struct distrOpt {
    int gszI, gszJ;         /* Global Size */
    int pnI, pnJ;          /* dim of 2d proc mesh */
    int nIter;             /* Number of iterations */
    int file_jump;         /* Number of iterations before new file gen*/
    int seed;
    int save_file;          /* save to file?*/
    int lIsz;
    int lJsz;
    char program_name[100];

};
struct neighbours {
    int left;
    int right;
    int down;
    int up;
};


void distr_borders(bool **board, neighbours nbr, MPI_Comm communicator, distrOpt &options) {


    static MPI_Datatype mpi_vect = MPI_DATATYPE_NULL;

    int lIsz = (options.gszI - 1) / options.pnI;
    int lJsz = (options.gszJ - 1) / options.pnJ;


    if (mpi_vect == MPI_DATATYPE_NULL) {
        MPI_Type_vector(lIsz, 1, lJsz + 2, MPI_CXX_BOOL, &mpi_vect);
        MPI_Type_commit(&mpi_vect);
    }

    MPI_Request req[4];
    MPI_Isend(&board[1][1], 1, mpi_vect, nbr.left, 0, communicator, req);
    MPI_Irecv(&board[1][0], 1, mpi_vect, nbr.left, 0, communicator, req + 2);
    MPI_Isend(&board[1][lJsz], 1, mpi_vect, nbr.right, 0, communicator, req + 1);
    MPI_Irecv(&board[1][lJsz + 1], 1, mpi_vect, nbr.right, 0, communicator, req + 3);

    MPI_Waitall(4, req, MPI_STATUSES_IGNORE);

    MPI_Isend(&board[1][0], lJsz + 2, MPI_CXX_BOOL, nbr.up, 0, communicator, req);
    MPI_Isend(&board[lIsz][0], lJsz + 2, MPI_CXX_BOOL, nbr.down, 0, communicator, req + 1);
    MPI_Irecv(&board[0][0], lJsz + 2, MPI_CXX_BOOL, nbr.up, 0, communicator, req + 2);
    MPI_Irecv(&board[lIsz + 1][0], lJsz + 2, MPI_CXX_BOOL, nbr.down, 0, communicator, req + 3);

    MPI_Waitall(4, req, MPI_STATUSES_IGNORE);

}


void initializeBoard(bool **board, distrOpt options) {
    int deadCellMultiplier = 2;
    srand(options.seed);

    for (int i = 0; i < options.lIsz; i++) {
        for (int j = 0; j < options.lJsz; j++) {
            board[i][j] = rand() % (deadCellMultiplier + 1) == 0;
        }
    }
}

void updateBoard(bool **board, distrOpt options) {
    const size_t rows = options.lIsz - 1;
    const size_t cols = options.lJsz - 1;
    std::vector<std::vector<int>> liveNeighbors(rows, std::vector<int>(cols, 0));

    //Count live neighbours
    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            if (board[i][j]) {
                for (int di = -1; di <= 1; ++di) {
                    for (int dj = -1; dj <= 1; ++dj) {
                        //Periodic boundary conditions
                        liveNeighbors[(i + di + rows) % rows][(j + dj + cols) % cols]++;
                    }
                }
                liveNeighbors[i][j]--; //Correction so that a cell does not consider itself as a live neighbor
            }
        }
    }

    //Update board
    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            board[i][j] = ((liveNeighbors[i][j] == 3) || (board[i][j] && liveNeighbors[i][j] == 2));
        }
    }
}

void writeBoardToFile(bool **board, size_t firstRow, size_t lastRow, size_t firstCol,
                      size_t lastCol, const std::string &fileName, int iteration, uint processID, distrOpt options) {
    //Open file
    std::ofstream outputFile(fileName + "_" + std::to_string(iteration) + "_" + std::to_string(processID) + ".gol");
    //Write metadata
    outputFile << std::to_string(firstRow) << " " << std::to_string(lastRow) << std::endl;
    outputFile << std::to_string(firstCol) << " " << std::to_string(lastCol) << std::endl;
    //Write data

    std::ostream_iterator<bool> outputIterator(outputFile, "\t");
    for (int i = 0; i < options.lIsz; i++) {
        std::vector<bool> line;
        for (int j = 1; j < options.lJsz; j++) {
            line.push_back(board[i][j]);
        }
        copy(line.begin(), line.end(), outputIterator);
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
    time(&rawtime);
    timeInfo = localtime(&rawtime);
    strftime(buffer, sizeof(buffer), "%Y-%m-%d-%H-%M-%S", timeInfo);
    std::string programName(buffer);

    //Generate main file
    std::ofstream outputFile(programName + ".gol");
    outputFile << std::to_string(rows) << " " << std::to_string(cols) << " " << std::to_string(iteration_gap) << " "
               << std::to_string(iterations) << " " << std::to_string(processes) << std::endl;
    outputFile.close();
    return programName;
}

int main(int argc, char *argv[]) {

    int rank, processors;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &processors);

    MPI_Datatype GOL_options;
    MPI_Aint displ[2] = {0, 10 * sizeof(int)};

    int lengths[2] = {10, 100};
    MPI_Datatype types[2] = {MPI_INT, MPI_CHAR};

    MPI_Type_create_struct(1, lengths, displ, types, &GOL_options);
    MPI_Type_commit(&GOL_options);

    distrOpt options;


    if (rank == 0) {
        if (argc != 5) {
            std::cout
                    << "This program should be called with four arguments! \nThese should be, the total number of rows; the total number of columns; the gap between saved iterations and the total number of iterations, in that order."
                    << std::endl;
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
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
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
            return 1;
        }
        float z = sqrt(processors);
        if (rows <= 0 || rows != cols || z != floor(z) || rows % (int) z != 0 || rows / (int) z < 4) {
            std::cout << "Illegal board size parameter combination!" << std::endl;
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
            return 1;
        }
        int mesh_dim = (int) z;
        std::string name = setUpProgram(rows, cols, iteration_gap, iterations, processors);
        std::cout << rank << ": Name: " << name << std::endl;
        int seed = rand();
        options = {
                cols, rows,            /* Global Size */
                mesh_dim, mesh_dim,          /* dim of 2d proc mesh */
                iterations,             /* Number of iterations */
                iteration_gap,         /* Number of iterations before new file gen*/
                seed,
                0,           /* save to file?*/
                cols / mesh_dim, rows / mesh_dim
        };
        name.copy(options.program_name, name.size());
        options.program_name[name.size()] = '\0';

//        int processes = 1, processID = 0;
//        size_t firstRow = 0, lastRow = rows - 1, firstCol = 0, lastCol = cols - 1;
//        std::string programName = setUpProgram(rows, cols, iteration_gap, iterations, processes);

        std::cout << rank << "Broadcasting options" << std::endl;

    }

    MPI_Bcast(&options, 1, GOL_options, 0, MPI_COMM_WORLD);

    MPI_Type_free(&GOL_options);

    std::string program_name = options.program_name;
    int pcoory = rank / (options.pnI);
    int pcoorx = rank % (options.pnJ );

    int p_up ,p_down,p_right,p_left;
    int dims[2]={options.pnI,options.pnJ};
    int periods[2]={0,0};

    MPI_Comm comm;
    MPI_Cart_create( MPI_COMM_WORLD, 2, dims, periods, 1, &comm );
    MPI_Cart_shift(comm,0,1,&p_up,&p_down);
    MPI_Cart_shift(comm,1,1,&p_left,&p_right);

    neighbours nbr = {p_left, p_right, p_down, p_up};

    int left_c = pcoorx * ((options.gszI ) / options.pnI);
    int right_c = (pcoorx + 1) * ((options.gszI) / options.pnI);
    int down_c = pcoory * ((options.gszJ) / options.pnJ);
    int up_c = (pcoory + 1) * ((options.gszJ) / options.pnJ);

    int lIsz = (options.gszI) / options.pnI;
    int lJsz = (options.gszJ) / options.pnJ;


    std::cout << rank << ": Neighbours determined" << p_left << " , " << p_right << " , " << p_down << " , " << p_up
              << std::endl;
    //Build board
    bool **board;
    board = (bool **) malloc(lIsz + 2 * sizeof(bool *));
    board[0] = (bool *) malloc((lIsz + 2) * (lJsz + 2) * sizeof(bool));
    std::cout << rank << ": Setting Borders" << std::endl;
    MPI_Barrier(MPI_COMM_WORLD);
    for (int i = 1; i < lIsz + 2; i++) {
        board[i] = board[i - 1] + lJsz + 2;
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0) {std::cout << rank << ": succeeded in setting pointers board" << std::endl;
        std::this_thread::sleep_for(std::chrono::milliseconds(1000));
    }
    MPI_Barrier(MPI_COMM_WORLD);
    initializeBoard(board, options);
    std::cout << rank << ": succeeded in initializing board" << std::endl;
    MPI_Barrier(MPI_COMM_WORLD);
    std::cout << rank << ": Board initialized" << std::endl;
    distr_borders(board, nbr, comm, options);
    MPI_Barrier(MPI_COMM_WORLD);
    std::cout << rank << ": Borders distributed" << std::endl;
    options.save_file=0;

    //Do iteration
    if (options.save_file) {
        writeBoardToFile(board, up_c, down_c, left_c, right_c, program_name, 0, rank, options);
    }
    for (int a = 1; a <= options.nIter; ++a) {
        updateBoard(board, options);
        std::cout << rank << ": update succeeded" << std::endl;
        MPI_Barrier(MPI_COMM_WORLD);
        if (options.save_file||a % options.file_jump == 0) {
            writeBoardToFile(board, up_c, down_c, left_c, right_c, program_name, a, rank, options);
        }
        distr_borders(board, nbr, comm, options);
        MPI_Barrier(MPI_COMM_WORLD);

    }

}



