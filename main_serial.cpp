//
// Created by Arthur Decloedt on 08/01/2020.
//

/**
 * Serial implementatin of the game of life for the course Parallel Computing.
 * Emil Loevbak (emil.loevbak@cs.kuleuven.be)
 * First implementation: November 2019
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>
#include <vector>
#include <ctime>

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
    const size_t rows = options.lIsz;
    const size_t cols = options.lJsz;
    std::vector<std::vector<int>> liveNeighbors(rows+2, std::vector<int>(cols+2, 0));

    //Count live neighbours
    for (size_t i = 1; i <= rows; ++i) {
        for (size_t j = 1; j <= cols; ++j) {
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
    for (size_t i = 1; i <= rows; ++i) {
        for (size_t j = 1; j <=cols; ++j) {
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
    if (argc != 5) {
        std::cout << "This program should be called with four arguments! \nThese should be, the total number of rows; the total number of columns; the gap between saved iterations and the total number of iterations, in that order." << std::endl;
        return 1;
    }
    int rows, cols;
    int iteration_gap, iterations;
    try{
        rows = atoi(argv[1]);
        cols = atoi(argv[2]);
        iteration_gap = atoi(argv[3]);
        iterations = atoi(argv[4]);
    } catch (std::exception const & exc) {
        std::cout << "One or more program arguments are invalid!" << std::endl;
        return 1;
    }
    int processes = 1, processID = 0;
    size_t firstRow = 0, lastRow = rows-1, firstCol = 0, lastCol = cols-1;
    std::string programName = setUpProgram(rows, cols, iteration_gap, iterations, processes);

    //Build board
    bool **board;
    board = (bool **) malloc((rows + 2) * sizeof(bool *));
    if (!board) {
        std::cout << "Couldn't allocate!" << std::endl;
        return 1;
    }

    board[0] = (bool *) malloc((rows + 2) * (cols + 2) * sizeof(bool));
    std::cout << "allocated board!" << std::endl;

    for (int i = 1; i < rows + 2; i++) {
        board[i] = board[i - 1] + cols + 2;
    }

    int seed = rand();
    distrOpt options = {
            cols, rows,            /* Global Size */
            1, 1,          /* dim of 2d proc mesh */
            iterations,             /* Number of iterations */
            iteration_gap,         /* Number of iterations before new file gen*/
            seed,
            0,           /* save to file?*/
            cols / 1, rows / 1
    };
    std::string program_name = options.program_name;
    int pcoory = 0 / (options.pnI);
    int pcoorx = 0 % (options.pnJ );

    int left_c = pcoorx * options.lIsz;
    int right_c = (pcoorx + 1) * options.lIsz;
    int down_c = pcoory * options.lJsz;
    int up_c = (pcoory + 1) * options.lJsz;

    initializeBoard(board,options);
    //Do iteration
    writeBoardToFile(board, down_c, up_c, left_c, right_c, programName, 0, 0, options);
    for (int i=1; i<=iterations; ++i) {
        updateBoard(board,options);
        std::cout << "update!" << std::endl;
        if (i%iteration_gap == 0) {
            writeBoardToFile(board, down_c, up_c, left_c, right_c, programName, i, 0, options);
        }
    }
}