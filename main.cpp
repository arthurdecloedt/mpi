
#include <iostream>
#include <fstream>
#include <iterator>
#include <vector>
#include <ctime>
#include <mpi.h>
#include <cmath>
#include <string>
#include <thread>
#include <chrono>


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

    int lIsz = options.lIsz;
    int lJsz = options.lIsz;


    if (mpi_vect == MPI_DATATYPE_NULL) {
        MPI_Type_vector(lIsz, 1, lJsz + 2, MPI_CXX_BOOL, &mpi_vect);
        MPI_Type_commit(&mpi_vect);
    }

    MPI_Request req[4];
    MPI_Isend(&board[1][1], 1, mpi_vect, nbr.right, 0, communicator, req);
    MPI_Irecv(&board[1][0], 1, mpi_vect, nbr.right, 0, communicator, req + 2);
    MPI_Isend(&board[1][lJsz], 1, mpi_vect, nbr.left, 0, communicator, req + 1);
    MPI_Irecv(&board[1][lJsz + 1], 1, mpi_vect, nbr.left, 0, communicator, req + 3);

    MPI_Waitall(4, req, MPI_STATUSES_IGNORE);

    MPI_Isend(&board[1][0], lJsz + 2, MPI_CXX_BOOL, nbr.up, 0, communicator, req);
    MPI_Isend(&board[lIsz][0], lJsz + 2, MPI_CXX_BOOL, nbr.down, 0, communicator, req + 1);
    MPI_Irecv(&board[0][0], lJsz + 2, MPI_CXX_BOOL, nbr.up, 0, communicator, req + 2);
    MPI_Irecv(&board[lIsz + 1][0], lJsz + 2, MPI_CXX_BOOL, nbr.down, 0, communicator, req + 3);

    MPI_Waitall(4, req, MPI_STATUSES_IGNORE);

}


void initializeBoard(bool **board, distrOpt options, int rank) {
    int deadCellMultiplier = 2;
    srand(rank);
    for (int i = 1; i < options.lIsz + 1; i++) {
        for (int j = 1; j < options.lJsz + 1; j++) {
            board[i][j] = rand() % (deadCellMultiplier + 1) == 0;
        }
    }

}

static bool inline next(bool **board, size_t i, size_t j) {
    int sum;

    sum = board[i - 1][j - 1] + board[i - 1][j] +
          board[i - 1][j + 1] + board[i][j - 1] +
          board[i][j + 1] + board[i + 1][j - 1] +
          board[i + 1][j] + board[i + 1][j + 1];

    if (sum < 2 || sum > 3) return false;
    else if (sum == 3) return true;
    else return board[i][j];
}


void updateBoard(bool **board, bool **nboard, distrOpt options) {
    const size_t rows = options.lIsz;
    const size_t cols = options.lJsz;

    //Update
    for (size_t i = 1; i <= rows; ++i) {
        for (size_t j = 1; j <= cols; ++j) {
            nboard[i][j] = next(board, i, j);
        }
    }
}


void writeBoardToFile(bool **board, size_t firstRow, size_t lastRow, size_t firstCol,
                      size_t lastCol, const std::string &fileName, int iteration, uint processID, distrOpt options) {
    //Open file

    std::ofstream outputFile("/data/leuven/332/vsc33219/gol/" + fileName + "_" + std::to_string(iteration) + "_" +
                             std::to_string(processID) + ".gol");
    //Write metadata
    outputFile << std::to_string(firstRow) << " " << std::to_string(lastRow) << std::endl;
    outputFile << std::to_string(firstCol) << " " << std::to_string(lastCol) << std::endl;
    //Write data

    std::ostream_iterator<bool> outputIterator(outputFile, "\t");
    for (int i = 1; i <= options.lIsz; i++) {
        std::vector<bool> line;
        for (int j = 1; j <= options.lJsz; j++) {
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
    std::this_thread::sleep_for(std::chrono::milliseconds(1000));
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    std::string time_file;
    int rank, processors;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &processors);

    MPI_Datatype GOL_options;
    MPI_Aint displ[2] = {0, 10 * sizeof(int)};

    int lengths[2] = {10, 100};
    MPI_Datatype types[2] = {MPI_INT, MPI_CHAR};
    MPI_Type_create_struct(2, lengths, displ, types, &GOL_options);
    MPI_Type_commit(&GOL_options);

    distrOpt options;
    //***C++11 Style:***

    int first = 0;

    if (rank == 0) {
        if (argc != 5 && argc != 6 && argc != 7) {
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
        std::string name = setUpProgram(rows, cols, iteration_gap, iterations, processors);
        if (argc > 5) {
            time_file = argv[5];
        } else {
            time_file = name;
        }
        if (argc > 6) {
            first = atoi(argv[6]);

        }
        name.copy(options.program_name, name.size());
        options.program_name[name.size()] = '\0';
        std::cout << rank << ": " << options.program_name << std::endl;


//        int processes = 1, processID = 0;
//        size_t firstRow = 0, lastRow = rows - 1, firstCol = 0, lastCol = cols - 1;
//        std::string programName = setUpProgram(rows, cols, iteration_gap, iterations, processes);

        std::cout << rank << "Broadcasting options" << std::endl;

    }
    MPI_Bcast(&options, 1, GOL_options, 0, MPI_COMM_WORLD);


    MPI_Type_free(&GOL_options);
    std::string name = options.program_name;

    int p_up, p_down, p_right, p_left;
    int dims[2] = {options.pnI, options.pnJ};

    MPI_Dims_create(processors, 2, dims);
    int periods[2] = {0, 0};

    MPI_Comm comm;
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 1, &comm);
    MPI_Cart_shift(comm, 0, 1, &p_up, &p_down);
    MPI_Cart_shift(comm, 1, 1, &p_left, &p_right);
    int coor[2];
    MPI_Cart_coords(comm, rank, 2, coor);
    int pcoory = coor[1];
    int pcoorx = coor[0];

    neighbours nbr = {p_left, p_right, p_down, p_up};
    int left_c = pcoory * options.lIsz;
    int right_c = (pcoory + 1) * options.lIsz - 1;
    int down_c = pcoorx * options.lJsz;
    int up_c = (pcoorx + 1) * options.lJsz - 1;

    int lIsz = (options.gszI) / options.pnI;
    int lJsz = (options.gszJ) / options.pnJ;

    std::cout << rank << ": Neighbours determined" << p_left << " , " << p_right << " , " << p_down << " , " << p_up
              << std::endl;
    //Build board
    bool **board, **board2, **tmp, *board_data, *board_data2;
    board = (bool **) malloc((lIsz + 2) * sizeof(bool *));
    board2 = (bool **) malloc((lIsz + 2) * sizeof(bool *));
    board_data = (bool *) malloc((lIsz + 2) * (lJsz + 2) * sizeof(bool));
    board_data2 = (bool *) malloc((lIsz + 2) * (lJsz + 2) * sizeof(bool));


    for (int i = 0; i < lIsz + 2; i++) {
        board[i] = &(board_data[i * (lJsz + 2)]);
        board2[i] = &(board_data2[i * (lJsz + 2)]);
    }

    initializeBoard(board, options, rank);
    std::cout << rank << ": succeeded in initializing board" << std::endl;
    MPI_Barrier(MPI_COMM_WORLD);
    std::cout << rank << ": Board initialized" << std::endl;
    distr_borders(board, nbr, comm, options);
    //Do iteration
//
//    writeBoardToFile(board, down_c, up_c, left_c, right_c, name, 0, rank, options);

    MPI_Barrier(MPI_COMM_WORLD);

    std::chrono::steady_clock::time_point check1 = std::chrono::steady_clock::now();

    for (int a = 1; a <= options.nIter; ++a) {
        updateBoard(board, board2, options);

        tmp = board;
        board = board2;
        board2 = tmp;
        MPI_Barrier(MPI_COMM_WORLD);

        distr_borders(board, nbr, comm, options);
        if (options.save_file && a % options.file_jump == 0) {
            writeBoardToFile(board, down_c, up_c, left_c, right_c, name, a, rank, options);
//            writeBoardToFile(board, left_c, right_c, down_c, up_c, name, a, rank, options);
        }

    }
    free(board_data);
    free(board);
    free(board_data2);
    free(board2);
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    MPI_Barrier(MPI_COMM_WORLD);
    long int local = std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count();
    long int local_nosetup = std::chrono::duration_cast<std::chrono::microseconds>(end - check1).count();
    long int local_setup = std::chrono::duration_cast<std::chrono::microseconds>(check1 - begin).count();
    long int global_sum;
    long int global_sum_nosetup;
    long int global_sum_setup;

    MPI_Reduce(&local, &global_sum, 1, MPI_LONG, MPI_SUM, 0,
               MPI_COMM_WORLD);
    MPI_Reduce(&local_nosetup, &global_sum_nosetup, 1, MPI_LONG, MPI_SUM, 0,
               MPI_COMM_WORLD);
    MPI_Reduce(&local_setup, &global_sum_setup, 1, MPI_LONG, MPI_SUM, 0,
               MPI_COMM_WORLD);

    if (rank == 0) {
        long int ttl_single = local;
        long int ttl_sum = global_sum;
        long int avg = global_sum / processors;
        long int avg_nosetup = global_sum_nosetup / processors;
        long int avgsetup = global_sum_setup / processors;

        std::ofstream outfile;
        outfile.open(time_file + "_detailed.out", std::ios_base::app); // append instead of overwrite
        outfile << "Timing results: milliseconds " << std::endl;

        outfile << "size:" << options.gszI << " by " << options.gszJ << std::endl;

        outfile << processors << " Processors" << std::endl;
        outfile << "Full (with setup) " << std::endl;
        outfile << "Single time (rank 0): " << ttl_single << "ms" << std::endl;
        outfile << "Avg single time: " << avg << "ms" << std::endl;
        outfile << "Summed time: " << ttl_sum << "ms" << std::endl;
        outfile << "Without setup " << std::endl;
        outfile << "Single time (rank 0): " << local_nosetup << "ms" << std::endl;
        outfile << "Avg single time: " << avg_nosetup << "ms" << std::endl;
        outfile << "Summed time: " << global_sum_nosetup << "ms" << std::endl;
        outfile << "Setup " << std::endl;
        outfile << "Single time (rank 0): " << local_setup << "ms" << std::endl;
        outfile << "Avg single time: " << avgsetup << "ms" << std::endl;
        outfile << "Summed time: " << global_sum_setup << "ms" << std::endl;
        outfile << "___________________________________________________ " << std::endl << std::endl;
        outfile.close();
        std::ofstream outfile2;
        outfile2.open(time_file + "_compact.csv", std::ios_base::app); // append instead of overwrite
        if (first != 0) {
            outfile2 << "X" << "," << "Y" << "," << "#P" << "," << "full single" << "," << "full avg" << ","
                     << "full sum" << "," << "nosetup single" << "," << "nosetup avg" << "," << "nosetup sum" << ","
                     << "setup single " << "," << "setup avg " << "," << "setup sum " << std::endl;
        }
        outfile2 << options.gszI << "," << options.gszJ << "," << processors << "," << ttl_single << "," << avg << ","
                 << ttl_sum << "," << local_nosetup << "," << avg_nosetup << "," << global_sum_nosetup << ","
                 << local_setup << "," << avgsetup << "," << global_sum_setup << std::endl;
        outfile2.close();
    }
    std::cout << rank << ": all succeeded" << std::endl;
    MPI_Finalize();
    return 0;
}



