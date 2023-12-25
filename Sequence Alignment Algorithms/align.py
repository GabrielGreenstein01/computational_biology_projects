
"""
STUDENT: Gabriel Greenstein
NETID: gbg222
STUDENT ID: 06605613

Collaborator(s): Upamanyu Dass-Vattam (udvattam)

Usage: python align.py input_file output_file

"""

# Import necessary packages
import sys
import numpy as np

#### ------ USEFUL FUNCTIONS ------- ####
def fuzzy_equals(a, b):
    """
    Checks if two floating point numbers are equivalent.
    """
    epsilon = 10**(-6) 
    return (abs(a - b) < epsilon)
    

#### ------- CLASSES ------- ####

class MatchMatrix(object):
    """
    Match matrix class stores the scores of matches in a data structure
    """
    def __init__(self):

        self.a = ""
        self.b = ""
        self.score = 0
        self.dictionary = {}

    def set_score(self, a, b, score):
        """
        Updates or adds a score for a specified match

        Inputs:
           a = the character from sequence A
           b = the character from sequence B
           score = the score to set it for
        """

        self.dictionary[(str(a),str(b))] = float(score)

    def get_score(self, a, b):
        """
        Returns the score for a particular match, where a is the
        character from sequence a and b is from sequence b.

        Inputs:
           a = the character from sequence A
           b = the character from sequence B
        Returns:
           the score of that match
        """

        return self.dictionary[(str(a),str(b))]


class Index(object):
    """
    Index stores the pointers and scores for a given index in the ScoreMatrix
    """
    def __init__(self):

        self.score = 0;
        self.pointers = []

    def add_pointer(self, name, row, col):
        """
        Adds a pointer to a given index;
        Contains the matrix to which it points and the index it points to (row,col)
        Input: matrix name, row, col to which the index points to (stored in the following form: [(matrix name, row, col),…………])
        """
        self.pointers.append(tuple((name,row,col)))

    def print_pointers(self):
        """
        Return: an index's pointers 
        """
        return self.pointers

    def index_score(self, scoreInput):
        """
        Sets an index's score
        Input: score for the given index
        """
        self.score = scoreInput

    def print_index_score(self):
        """
        Return: an index's score
        """
        return self.score


class ScoreMatrix(object):
    """
    Object to store a score matrix, which generated during the alignment process. The score matrix consists of a 2-D array of
    ScoreEntries that are updated during alignment and used to output the maximum alignment.
    """

    def __init__(self, name, nrow, ncol):

        self.name = str(name)
        self.nrow = int(nrow)
        self.ncol = int(ncol)
        self.score_matrix = []

        # Initialize matrix with Index objects
        for i in range(0,nrow):
            row = []
            for j in range(0,ncol):
                row.append(Index())

            self.score_matrix.append(row)

    def get_score(self, row, col):
        """
        Return: score of a matrix's index; uses Index object and print_index_score()
        """
        return self.score_matrix[row][col].print_index_score()
   
    def set_score(self, row, col, score):    
        """
        Sets the score of a matrix's index; uses Index object and index_score()
        """
        self.score_matrix[row][col].index_score(float(score))

    def get_pointers(self, row, col):
        """
        Returns the indices of the entries that are pointed to
        This should be formatted as a list of tuples:
         ex. [(matrix name, row, col),…………]
        """
        return self.score_matrix[row][col].print_pointers()

    def set_pointers(self, curRow, curCol, pointMatrix, pointRow, pointCol):
        """
        Adds a pointer to a given index in a matrix; uses Index object and add_pointer()
        """
        self.score_matrix[curRow][curCol].add_pointer(pointMatrix,pointRow,pointCol)

    def print_scores(self):
        """
        Returns a nicely formatted string containing the scores in the score matrix. Use this for debugging!

        Example:
        M=
            0.0, 0.0, 0.0, 0.0, 0.0
            0.0, 1.0, 0.0, 0.0, 0.0
            0.0, 1.0, 1.0, 1.0, 1.0
            0.0, 0.0, 1.0, 1.0, 1.0
            0.0, 0.0, 2.0, 2.0, 1.0
            0.0, 0.0, 1.0, 2.0, 3.0

        """

        # Create matrix to store just the scores of each index so we can display it nicely
        matrixPrint = []

        # Get the score of every index in a given matrix and append to matrixPrint()
        for i in range(0,self.nrow):
            row = []
            for j in range(0,self.ncol):
                row.append(self.score_matrix[i][j].print_index_score())
            matrixPrint.append(row)

        print(np.matrix(matrixPrint))

        return

    def print_pointers(self):
        """
        Returns a nicely formatted string containing the pointers for each entry in the score matrix. Use this for debugging!
        """
        pointerMatrix = []

        # Get list of pointers for each index in a given matrix
        for i in range(0,self.nrow):
            for j in range(0,self.ncol):
                row = []
                row.append(str("Index: (" + str(i) + ", " + str(j) + ")"))
                row.append(str(self.score_matrix[i][j].print_pointers()))
                pointerMatrix.append(row)

        for row in pointerMatrix:
            print(*row)

        return


class AlignmentParameters(object):
    """
    Object to hold a set of alignment parameters from an input file.
    """

    def __init__(self):
        # default values for variables that are filled in by reading
        # the input alignment file
        self.seq_a = ""
        self.seq_b = ""
        self.global_alignment = False 
        self.dx = 0
        self.ex = 0
        self.dy = 0
        self.ey = 0
        self.alphabet_a = "" 
        self.alphabet_b = ""
        self.len_alphabet_a = 0
        self.len_alphabet_b = 0
        self.match_matrix = MatchMatrix()

    def load_params_from_file(self, input_file): 
        """
        Reads the parameters from an input file and stores in the object

        Input:
           input_file = specially formatted alignment input file
        """

        f = open(input_file, "r")
        contents = f.readlines()

        # Remove '\n' from each line
        for i in range(0,len(contents)):
            line = contents[i]
            contents[i] = line.strip("\n")

        # Define parameters based on input file
        self.seq_a = contents[0]
        self.seq_b = contents[1]

        # Update alignment type
        if contents[2] == '0':
            self.global_alignment = True
        else:
            self.global_alignment = False

        # Obtain dx, ex, dy, ey as array
        gap_param = contents[3].split(" ")

        # Update dx, ex, dy, ey
        self.dx = float(gap_param[0])
        self.ex = float(gap_param[1])
        self.dy = float(gap_param[2])
        self.ey = float(gap_param[3])


        # Update len_alphabet_a, alphabet_a, len_alphabet_b, alphabet_b
        self.len_alphabet_a = int(contents[4])
        self.alphabet_a = contents[5]
        self.len_alphabet_b = int(contents[6])
        self.alphabet_b = contents[7]

        # Remove first eight rows so that what remains are the match matrix lines
        contents = contents[8:]

        # Iterate through each line of match matrix and update self.match_matrix accordingly
        for row in contents:
            if len(row) != 0:
                row = row.split(" ")
                # row = row.split("\t") # Use for quiz_input/quiz5_pam40_human_chimp.input

                row = list(filter(None, row))
                self.match_matrix.set_score(str(row[2]),str(row[3]),float(row[4]))

        f.close()


class Align(object):
    """
    Object to hold and run an alignment; running is accomplished by using "align()"
    """

    def __init__(self, input_file, output_file):
        """
        Input:
            input_file = file with the input for running an alignment
            output_file = file to write the output alignments to
        """
        self.input_file = input_file
        self.output_file = output_file
        self.align_params = AlignmentParameters()

        # Initialize to ScoreMatrix of size (0,0); will update in populate_score_matrices()
        self.m_matrix = ScoreMatrix("m",0,0)
        self.ix_matrix = ScoreMatrix("ix",0,0)
        self.iy_matrix = ScoreMatrix("iy",0,0)

    def align(self):
        """
        Main method for running alignment.
        """
        # load the alignment parameters into the align_params object
        self.align_params.load_params_from_file(self.input_file)

        # Populate the score matrices based on the input parameters
        self.populate_score_matrices()

        # Perform a traceback on the possible alignments (i.e., paths)
        [score, paths] = self.traceback_helper()

        # Get best alignment score
        score = round(score,1)

        # Call write_output() to write .output file
        self.write_output(score, paths)

    def populate_score_matrices(self):
        """
        Method to populate the score matrices based on the data in align_params.
        Should call update(i,j) for each entry in the score matrices
        """

        # Initialize m_matrix, ix_matrix, iy_matrix to have dimensions (len(seq_a) + 1, len(seq_b) + 1)
        # The +1 is important because that's how the algorithm sets up the matrices (reference diagram)
        self.m_matrix = ScoreMatrix("M", len(self.align_params.seq_a) + 1, len(self.align_params.seq_b) + 1)
        self.ix_matrix = ScoreMatrix("Ix", len(self.align_params.seq_a) + 1, len(self.align_params.seq_b) + 1)
        self.iy_matrix = ScoreMatrix("Iy", len(self.align_params.seq_a) + 1, len(self.align_params.seq_b) + 1)

        # Call update() on each index aside from first row and first column
        for i in range(1, len(self.align_params.seq_a) + 1):
            for j in range(1, len(self.align_params.seq_b) + 1):
                self.update(i,j)

    def update(self, row, col):
        """
        Method to update the matrices at a given row and column index.

        Input:
           row = the row index to update
           col = the column index to update
        """
        # Update each matrix at the given index
        self.update_m(row, col)
        self.update_ix(row, col)
        self.update_iy(row, col)

    def update_m(self, row, col):
        """
        Function to update the scores of matrix M and add pointers to each index
        Input: row, col in score matrix m to be updated 
        """
        i = row
        j = col

        # Find characters to compare in match_matrix
        a = self.align_params.seq_a[i - 1]
        b = self.align_params.seq_b[j - 1]

        # Get score of a,b
        match_score = self.align_params.match_matrix.get_score(a,b)

        # Compute 3 possible scores for M(i,j)
        m = self.m_matrix.get_score(i-1,j-1) + match_score
        ix = self.ix_matrix.get_score(i-1,j-1) + match_score
        iy = self.iy_matrix.get_score(i-1,j-1) + match_score

        # Define nums and add the 3 options to it
        nums = []
        nums.extend([round(m,6),round(ix,6),round(iy,6)])

        # Get the max value of m, ix, iy 
        max_val = np.amax(nums)

        # Find which matrices (of m, ix, iy) the max value was generate from
        max_indices = np.argwhere(nums == max_val)
        max_indices = max_indices.flatten().tolist()

        # LOCAL ALIGNMENT ONLY: if the max value < 0, set max_val = 0
        # Recall that there are no negative numbers in any score matrix for local alignment
        if self.align_params.global_alignment == False and max_val < 0:
            max_val = 0

        # Set the score for m(i,j)
        self.m_matrix.set_score(i,j,max_val)

        # Add pointers for that new index based on which matrix and index yielded the max value
        for val in range(0,len(max_indices)):

            if max_indices[val] == 0:
                self.m_matrix.set_pointers(i,j,"m",i-1,j-1)
            if max_indices[val] == 1:
                self.m_matrix.set_pointers(i,j,"ix",i-1,j-1)
            if max_indices[val] == 2:
                self.m_matrix.set_pointers(i,j,"iy",i-1,j-1)

    def update_ix(self, row, col):
        """
        Function to update the scores of matrix Ix and add pointers to each index
        Input: row, col in score matrix Ix to be updated 
        """
        i = row
        j = col

        # Compute 2 possible scores for Ix(i,j)
        m = self.m_matrix.get_score(i-1,j) - self.align_params.dy
        ix = self.ix_matrix.get_score(i-1,j) - self.align_params.ey

        ##### Follows the same process as in update_m()
        nums = []
        nums.extend([round(m,6),round(ix,6)])

        max_val = np.amax(nums)

        max_indices = np.argwhere(nums == max_val)
        max_indices = max_indices.flatten().tolist()

        if self.align_params.global_alignment == False and max_val < 0:
            max_val = 0

        self.ix_matrix.set_score(i,j,max_val)

        for val in range(0,len(max_indices)):

            if max_indices[val] == 0:
                self.ix_matrix.set_pointers(i,j,"m",i-1,j)
            if max_indices[val] == 1:
                self.ix_matrix.set_pointers(i,j,"ix",i-1,j)

        return

    def update_iy(self, row, col):
        """
        Function to update the scores of matrix Iy and add pointers to each index
        Input: row, col in score matrix Iy to be updated 
        """
        i = row
        j = col

        # Compute 2 possible scores for Iy(i,j)
        m = self.m_matrix.get_score(i,j-1) - self.align_params.dx
        iy = self.iy_matrix.get_score(i,j-1) - self.align_params.ex

        ##### Follows the same process as in update_m()
        nums = []
        nums.extend([round(m,6),round(iy,6)])

        max_val = np.amax(nums)

        max_indices = np.argwhere(nums == max_val)
        max_indices = max_indices.flatten().tolist()

        if self.align_params.global_alignment == False and max_val < 0:
            max_val = 0

        self.iy_matrix.set_score(i,j,max_val)

        for val in range(0,len(max_indices)):

            if max_indices[val] == 0:
                self.iy_matrix.set_pointers(i,j,"m",i,j-1)
            if max_indices[val] == 1:
                self.iy_matrix.set_pointers(i,j,"iy",i,j-1)

        return

    def find_traceback_start(self):
        """
        Finds the location to start the traceback.

        Returns:
            (max_val, [(row,col),…………]) = Tuple[float, Set[Tuple[int, int]]]
        """
        # Find starting index/indices for LOCAL ALIGNMENT
        if self.align_params.global_alignment == False:
            
            # Initialize max_val to be tuple of the form (max_val, [(row, column)])
            max_val = tuple((0,set(tuple((0,0)))))

            # Search through all indices of m_matrix
            for i in range(0,self.m_matrix.nrow):
                for j in range(0,self.m_matrix.ncol):

                    # Get score of index
                    score = self.m_matrix.get_score(i,j)

                    # If the score is equal to max score found, add tuple to the set in max_val (i.e., max_val[1])
                    if score == round(max_val[0],6):
                        max_val[1].add(tuple((i,j)))

                    # If the score is > than max score found, reassign max_val to a new tuple containing the max score and index
                    if score > round(max_val[0],6):
                        max_val = tuple((score,set()))
                        max_val[1].add(tuple((i,j)))

            return max_val
                        
        else: # same comments as if statement

            # Initialize max_val to be tuple of the form (max_val, [(row, column)])
            max_val = tuple((0,set(tuple((0,0)))))

            # Search through last row and column of m_matrix:
            # Search through last column INCLUDING last row
            for i in range(0,self.m_matrix.nrow):

                score = self.m_matrix.get_score(i,self.m_matrix.ncol - 1)

                if score == round(max_val[0],6):
                    max_val[1].add(tuple((i,self.m_matrix.ncol - 1)))

                if score > round(max_val[0],6):
                    max_val = tuple((score,set()))
                    max_val[1].add(tuple((i,self.m_matrix.ncol - 1)))

            # Search through last row UP TO last column
            # (don't look at last column to avoid looking at bottom-right index again)
            for j in range(0,self.m_matrix.ncol - 1):

                score = self.m_matrix.get_score(self.m_matrix.nrow - 1,j)

                if score == round(max_val[0],6):
                    max_val[1].add(tuple((self.m_matrix.nrow - 1,j)))

                if score > round(max_val[0],6):
                    max_val = tuple((score,set()))
                    max_val[1].add(tuple((self.m_matrix.nrow - 1,j)))

            return max_val



    def traceback_helper(self):
        """
        Helper function to call our traceback recursively and keep track of all possible alignments starting from the max score
        Return: the max score (retrieved from start_traceback()) and the paths (i.e., alignments) originating from the starting points
        """
        # Get list of locations to start traceback, which had the max score for the given alignment criteria
        start_traceback = self.find_traceback_start()

        # Define global variable to keep track of all alignments, i.e. paths
        global paths
        paths = []

        # Get max score from find_traceback_start()
        max_val = start_traceback[0]

        # Define a list which will contain only the indices in start_traceback
        starting_points = []

        # Update start_points with the indices in start_traceback
        for point in start_traceback[1]:
            starting_points.append(tuple(("m",point[0],point[1])))

        # Call traceback()
        self.traceback([], starting_points)

        return [max_val, paths]

    def traceback(self, path, pointers): 
        """
        This is where we traceback a starting pointer to retrieve all possible alignments.
        Input: current path, pointers to add to new path
        Further, we have path as an input so we can recursively add new pointers to it and terminate it when the alignment is done
        """
        # Iterate through all pointers
        for i in range(0,len(pointers)):

            # Get the matrix name, col, and row for each pointer
            # Each pointer in the list has the form (max_val, max_matrix, row, column)
            row = pointers[i][1]
            col = pointers[i][2]
            matrix_name = pointers[i][0]

            # Base case(s): If pointer points to an index with no pointers (i.e., any index in the first row/col)
            if matrix_name == "m" and len(self.m_matrix.get_pointers(row,col)) == 0:
         
                paths.append(path)

                return

            elif matrix_name == "ix" and len(self.ix_matrix.get_pointers(row,col)) == 0:

                paths.append(path)

                return

            elif matrix_name == "iy" and len(self.iy_matrix.get_pointers(row,col)) == 0:

                paths.append(path)

                return

            else: # If pointer points to something else

                # Iterate through all possible matrices
                if matrix_name == "m":

                    # LOCAL ALIGNMENT: terminates at first 0 it encounters in the given score matrix 
                    if self.align_params.global_alignment == False and self.m_matrix.get_score(row,col) == 0:

                        paths.append(path)

                        return

                    # Otherwise, get the pointers
                    m_pointers = self.m_matrix.get_pointers(row,col)

                    # Create a new path (MUST use .copy() to allocate new memory) and add pointer
                    new_path = path.copy()
                    new_path.append(tuple((matrix_name,row,col)))

                    # Recursively call on traceback using the updated path and the new pointers
                    self.traceback(new_path,list((m_pointers)))

                    # delete newly allocated memory
                    del new_path

                if matrix_name == "ix": # Same comments as m

                    if self.align_params.global_alignment == False and self.ix_matrix.get_score(row,col) == 0:

                        paths.append(path)

                        return

                    ix_pointers = self.ix_matrix.get_pointers(row,col)

                    new_path = path.copy()
                    new_path.append(tuple((matrix_name,row,col)))

                    self.traceback(new_path,list((ix_pointers)))

                    del new_path

                if matrix_name == "iy": # Same comments as m

                    if self.align_params.global_alignment == False and self.iy_matrix.get_score(row,col) == 0:

                        paths.append(path)

                        return

                    iy_pointers = self.iy_matrix.get_pointers(row,col)

                    new_path = path.copy()
                    new_path.append(tuple((matrix_name,row,col)))

                    self.traceback(new_path,list((iy_pointers)))

                    del new_path

    def write_output(self, score, paths):
        """
        This writes the output file in the desired format
        Input: max score for which we start our traceback from (abbreviated 'score'), and the paths (i.e., alignments)
        originating from the traceback start
        """
        # Create new file whose name is the output filename provided in the command line
        f = open(self.output_file, "w")

        # Write the score of the alignment
        f.write(str(score) + "\n")

        # Iterate through each of the alignments
        for i in range(0,len(paths)):

            # Define empty aligned sequence which we will populate
            seq_a = ""
            seq_b = ""

            path = paths[i]

            # Iterate through each tuple ((matrix_name, row, col)) in the path
            for j in range(0,len(path)):

                # Remove the last element from path
                elem = path.pop()

                # Get the matrix name, row, and column
                matrix_name = elem[0]
                row = elem[1]
                col = elem[2]

                # Update the strings according to what M, Ix, and Iy do
                if matrix_name == "m":
                    seq_a = seq_a + self.align_params.seq_a[row - 1]
                    seq_b = seq_b + self.align_params.seq_b[col - 1]

                if matrix_name == "ix":
                    seq_a =  seq_a + self.align_params.seq_a[row - 1]
                    seq_b = seq_b + "_" 

                if matrix_name == "iy":
                    seq_a = seq_a + "_"
                    seq_b = seq_b + self.align_params.seq_b[col - 1]

            # Write sequences to new document
            f.write("\n" + seq_a + "\n")
            f.write(seq_b + "\n")

        f.close()


def main():
    """
    Function is called when the file is called in the terminal
    """
    # check that the file is being properly used
    if (len(sys.argv) !=3):
        print("Please specify an input file and an output file as args.")
        return
        
    # input variables
    input_file = sys.argv[1]
    output_file = sys.argv[2]

    # create an align object and run
    align = Align(input_file, output_file)
    align.align()

if __name__=="__main__":
    main()
