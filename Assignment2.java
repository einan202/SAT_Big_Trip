
public class Assignment2 {

	// Task 1
	public static boolean isSquareMatrix(boolean[][] matrix) {

		boolean isSquare = notNull(matrix);                 // no length or cells checks if null array

		for (int x = 0; isSquare && (x < matrix.length); x = x + 1) {

			int inArrayLength = matrix[x].length;
			isSquare = (inArrayLength == matrix.length);    // each cell's length compare to main's length

		}

		return isSquare;
	}


	public static boolean notNull(boolean[][] matrix) {
		// insures main array and its inner arrays not null

		boolean notNull = (matrix != null);

		for (int x = 0; notNull && (x < matrix.length); x = x + 1) {

			notNull = (matrix[x] != null);

		}

		return notNull;
	}


	// Task 2
	public static boolean isSymmetricMatrix(boolean[][] matrix) {

		boolean isSymmetric = notNull(matrix);                                       // no cells check if null array

		for (int i = 0; isSymmetric && (i < matrix.length - 1); i = i + 1) {         // keeps 'j in bounds

			for (int j = i + 1; isSymmetric & (j < matrix.length); j = j + 1) {     // no double checks for same partners

				isSymmetric = (matrix[i][j] == matrix[j][i]);

			}

		}

		return isSymmetric;
	}


	// Task 3
	public static boolean isAntiReflexiveMatrix(boolean[][] matrix) {

		boolean isAntiReflexive = notNull(matrix);        // no cells check if null array

		for (int x = 0; isAntiReflexive && (x < matrix.length); x = x + 1) {

			isAntiReflexive = (matrix[x][x] == false);    // negative value means no flight

		}

		return isAntiReflexive;
	}


	// Task 4
	public static boolean isLegalInstance(boolean[][] matrix) {

		boolean isLegal = notNull(matrix);                // null matrix not legal - no other checks

		if (isLegal) {

			isLegal =
					isSquareMatrix(matrix) &
							isSymmetricMatrix(matrix) &
							isAntiReflexiveMatrix(matrix);

		}

		return isLegal;
	}


	// Task 5
	public static boolean isPermutation(int[] array) {

		boolean isPermutation = true;                           // assumes array is'nt null

		for (int i = 0; isPermutation & (i < array.length); i = i + 1) {

			isPermutation =
					valueInBounds(array, i) &&
							IndexAppears(array, i) &&           // no out-of-bounds values exist
							countAppears(array, i) == 1;        // each city appears once

		}

		return isPermutation;
	}


	public static boolean valueInBounds(int[] array, int i) {
		// checks each cell's value is in bounds

		return !(array[i] < 0 | array[i] > array.length - 1);
	}


	public static boolean IndexAppears(int[] array, int i) {
		// insures every city appears in the tour

		boolean found = false;

		for (int j = 0; !found & (j < array.length); j = j + 1) {

			found = (array[j] == i);

		}

		return found;
	}


	public static int countAppears(int[] array, int i) {
		// counts the times an index appears in tour

		int counter = 0;

		for (int j = 0; j < array.length; j = j + 1) {

			if (array[j] == i) counter = counter + 1;

		}

		return counter;
	}


	// Task 6
	public static boolean hasLegalSteps(boolean[][] flights, int[] tour) {

		return tourExist(flights, tour);
	}


	public static boolean tourExist(boolean[][] flights, int[] tour) {
		// checks only relevant flights in tour

		boolean tourExist = true;

		for (int i = 0; i < tour.length & tourExist; i = i + 1) {

			int origin = tour[i];
			int destination = tourIndex(flights, tour, i);         // so last city destination is first

			tourExist =
					flightExist(flights, origin, destination);

		}

		return tourExist;
	}


	public static boolean flightExist(boolean[][] flights, int origin, int destination) {
		// checks flight exist between origin and destination cities

		boolean flightExist = true;

		for (int i = 0; i < flights[origin].length & flightExist; i = i + 1) {

			flightExist = (flights[origin][destination]); // true means flight exists

		}

		return flightExist;
	}


	public static int tourIndex(boolean[][] flights, int[] tour, int i) {
		// checks which destination is relevant one

		int destination;

		if (i != tour.length - 1)                       // different than the last
			destination = tour[i + 1];

		else {
			destination = tour[0]; 					    // last city back to first one
		}

		return destination;
	}


	// Task 7
	public static boolean isSolution(boolean[][] flights, int[] tour) {

		int n = flights.length;                        // assumes flights array is'nt null
		isSolutionException(tour, n);

		return (tour[0] == 0) &&                       // shortest to the longest check
				isPermutation(tour) &&
				hasLegalSteps(flights, tour);
	}


	public static void isSolutionException(int[] tour, int n) {
		// throws an exception for an illegal tour

		if (tour == null || tour.length != n)          // checks nullity before checking length
			throw new IllegalArgumentException("tour is null or it's length is'nt legal");

	}


	// Task 8
	public static int[][] atLeastOne(int[] vars) {

		int[][] cnf = new int[1][vars.length];         // 'or' between all vars is only one clause
		return atLeastFormula(vars, cnf);
	}


	public static int[][] atLeastFormula(int[] vars, int[][] cnf) {
		// creates the formula for minimum one true

		for (int x = 0; x < vars.length; x = x + 1) {

			cnf[0][x] = vars[x];                       // each cell in first (and only) array is a var

		}

		return cnf;
	}


	// Task 9
	public static int[][] atMostOne(int[] vars) {

		int combinations = ((vars.length) * (vars.length - 1)) / 2;    // couples possible combinations
		int[][] cnf = new int[combinations][2];                        // each condition is a clause

		return atMostFormula(vars, cnf, 0);                    // separate counter for clause index
	}


	public static int[][] atMostFormula(int[] vars, int[][] cnf, int counter) {
		// creates the formula for maximum one true

		for (int x = 0; x < vars.length; x = x + 1) {

			for (int y = x + 1; y < vars.length; y = y + 1) {          // prevents double combinations

				int[] clause = {-(vars[x]), -(vars[y])};               // party couples formula
				cnf[counter] = clause;
				counter = counter + 1;

			}

		}

		return cnf;
	}


	// Task 10
	public static int[][] exactlyOne(int[] vars) {

		int combinations = (vars.length * (vars.length - 1)) / 2;      // couples possible combinations
		int[][] cnf = new int[1 + combinations][];                     // conjunction of atLeast and atMost clauses

		return exactlyOneFormula(vars, cnf, combinations);
	}


	public static int[][] exactlyOneFormula(int[] vars, int[][] cnf, int combinations) {
		// creates one formula of atLeast and atMost together

		cnf[0] = atLeastOne(vars)[0];                       // insures minimum one true

		for (int x = 1; x <= combinations; x = x + 1) {

			cnf[x] = atMostOne(vars)[x - 1];                // contains all maximum one clauses

		}

		return cnf;
	}


	// Task 11
	public static boolean[] solveExactlyOneForEachSet(int[][] varSets) {

		SATSolver.init(varsMaxValue(varSets));               // enables SAT working with all set's vars

		for (int s = 0; s < varSets.length; s = s + 1)
			SATSolver.addClauses(exactlyOne(varSets[s]));    // each Formula for a set

		timeOutException(SATSolver.getSolution());           // exception for null solution
		return SATSolver.getSolution();
	}


	public static int maxValue(int[] varSet) {
		// finds max variable of one set

		int maxValue = varSet[0];                            // same as MinIndex but of max

		for (int x = 1; x < varSet.length; x = x + 1) {

			int currValue = varSet[x];
			if (currValue > maxValue) maxValue = currValue;

		}

		return maxValue;
	}


	public static int varsMaxValue(int[][] varSets) {
		// calculates max variable of all sets

		int maxValue = 0;
		for (int s = 1; s < varSets.length - 1; s = s + 1) {

			maxValue = maxValue(varSets[s]);            // max of a set

			if (maxValue(varSets[s + 1]) > maxValue)    // max of all sets
				maxValue = maxValue(varSets[s + 1]);

		}

		return maxValue;
	}


	public static void timeOutException(boolean[] solution) {
		// throws an exception for too long time check

		if (solution == null)
			throw new IllegalArgumentException("TIMEOUT");
	}


	// Task 12
	public static int[][] createVarsMap(int n) {

		int[][] map = new int[n][n];

		for (int i = 0; i < n; i = i + 1) {
			for (int j = 0; j < n; j = j + 1) {

				map[i][j] = (i * n) + j + 1;     // not repeating same vars values

			}
		}

		return map;
	}


	// Task 13
	public static int[][] oneCityInEachStep(int[][] map) {

		int steps = map.length;                                        // map size is n * n
		int oneStepFormula = 1 + ((steps * (steps - 1)) / 2);          // exactlyOne formula's length

		int[][] cnf = new int[steps * oneStepFormula][];               // length for all steps formula

		return exactlyMixOne(cnf, map, oneStepFormula);
	}


	public static int[][] exactlyMixOne(int[][] cnf, int[][] map, int oneStepFormula) {
		// creates one formula of exactlyOne each step

		int fullFormulaIndex = 0;

		for (int i = 0; i < map.length; i = i + 1) {

			for (int clause = 0; clause < oneStepFormula; clause = clause + 1) {

				cnf[fullFormulaIndex] = exactlyOne(map[i])[clause];        // map [i] is one step. for example: {1,2,3}
				fullFormulaIndex = fullFormulaIndex + 1;                   // separate formula index

			}

		}

		return cnf;
	}


	// Task 14
	public static int[][] fixSourceCity(int[][] map) {

		int[][] oneClause = {{map[0][0]}};                           // first var in first step
		return oneClause;
	}


	// Task 15
	public static int[][] eachCityIsVisitedOnce(int[][] map) {

		int steps = map.length;                                        // each step - one var
		int oneStepFormula = 1 + ((steps * (steps - 1)) / 2);          // exactly one's length

		int[][] cnf = new int[steps * oneStepFormula][];
		int[][] upsideMap = upsideMap(map);                            // flips map - now formula is between same cities

		return exactlyMixOne(cnf, upsideMap, oneStepFormula);
	}


	public static int[][] upsideMap(int[][] map) {
		// creates a flipped map

		int[][] upsidemap = new int[map.length][map.length];      // separate to not change original

		for (int i = 0; i < map.length; i = i + 1) {

			for (int j = 0; j < map[i].length; j = j + 1) {

				upsidemap[j][i] = (i * map.length) + j + 1;       // numerates vertically

			}
		}

		return upsidemap;
	}


	// Task 16
	public static int[][] noIllegalSteps(boolean[][] flights, int[][] map) {

		int [][] nonFlights = nonFlightsArray(flights);    // non-flight cities list
		int nonCouples = nonFlights.length;				   // non flights counter

		int [][] cnf = new int [(1 + 2 * (map.length - 1)) * nonCouples][];

		return createIllegalCNF(cnf, nonCouples, nonFlights, map);
	}


	public static int [][] createIllegalCNF(int [][] cnf, int nonCouples, int [][] nonFlights, int [][] map){
		// creates final noIllegalSteps formula

		int cnfIndx = 0;
		boolean cond = cnfIndx < cnf.length;         // in case of no non-flights

		for (int x = 0; cond && x < nonCouples; x = x + 1) {
			cond = cnfIndx < cnf.length;

			int i = 0, j = nonFlights[x][0], k = nonFlights[x][1];

			while (cond && i < map.length - 1) {

				createCouplesCNF(cnf, cnfIndx, map, i, j, k);
				cnfIndx = cnfIndx + 2;               // two options which first
				i = i + 1;
			}

			createLastsCNF(cnf, cnfIndx, map, i,k);  // back to origin
			cnfIndx = cnfIndx + 1;
		}

		return cnf;
	}


	public  static void createCouplesCNF(int [][] cnf, int cnfIndx ,int [][] map, int i, int j, int k) {
		// adds non-flights couples formula without last

		int c1 = map[i][j], c2 = map[i + 1][k];  // in case j is first
		int[] clauseA = {-c1, -c2};

		int c3 = map[i][k], c4 = map[i + 1][j];  // in case k is first
		int[] clauseB = {-c3, -c4};

		cnf[cnfIndx] = clauseA;
		cnf[cnfIndx + 1] = clauseB;
	}


	public static void createLastsCNF(int [][] cnf, int cnfIndx, int [][] map, int i, int k){
		// adds back to origin city clause

		if (i == map.length - 1) { 	// last in tour

			int [] clauseC = {-map[i][k]};
			cnf[cnfIndx]  = clauseC;

		}
	}


	public static int [][] nonFlightsArray(boolean[][] flights) {
		// returns an array of cities with no flight

		int nonFlights = nonFlights(flights);
		int [][] arr = new int[nonFlights][2];  // flight between two only

		int counter = 0;
		boolean c1 = counter < nonFlights;

		for (int j = 0; c1 && j < flights.length - 1; j = j + 1) {

			c1 = counter < nonFlights;		    // for no out of bound
			for (int k = j + 1; c1 && k < flights.length; k = k + 1) {

				if (j != k & !flightExist(flights, j, k)){

					arr[counter][0] = j;	    // first city
					arr[counter][1] = k;	    // second city
					counter = counter + 1;
				}
			}
		}

		return arr;
	}


	public static int nonFlights(boolean[][] flights) {
		// returns number of couples with no flight

		int nonFlights = 0;

		for (int j = 0; j < flights.length - 1; j = j + 1) {
			for (int k = j + 1; k < flights.length; k = k + 1) {

				if (j != k & !flightExist(flights, j, k)) {  // assumes Anti-reflexive

					nonFlights = nonFlights + 1;

				}
			}
		}

		return nonFlights;
	}


	// Task 17
	public static void encode(boolean[][] flights, int[][] map) {

		if (flights == null || !isLegalInstance(flights))                         // checks null before calling function
			throw new IllegalArgumentException("flights array is illegal");

		if (map == null || !validMap(map, createVarsMap(flights.length)))     // assumes flights check passed ok
			throw new IllegalArgumentException("map does not match");

		insertClauses(flights, map);
	}


	public static boolean validMap(int[][] map, int[][] legalMap) {
		// compares a valid map of flights with input one

		boolean validMap = (map.length == legalMap.length);

		for (int x = 0; validMap & x < map.length; x = x + 1) {
			validMap = (map[x] != null && map[x].length == legalMap[x].length); // checks null before lengths

			for (int y = 0; validMap & y < map.length; y = y + 1) {
				validMap = (map[x][y] == legalMap[x][y]);                        // same var number
			}
		}

		return validMap;
	}


	public static void insertClauses(boolean[][] flights, int[][] map) {
		// already initialized, adding relevant clauses to solver

		SATSolver.addClauses(oneCityInEachStep(map));
		SATSolver.addClauses(eachCityIsVisitedOnce(map));
		SATSolver.addClauses(fixSourceCity(map));
		SATSolver.addClauses(noIllegalSteps(flights, map));
	}


	// Task 18
	public static int[] decode(boolean[] assignment, int[][] map) {

		int n = map.length;                      // assumes map is legal

		if (assignment.length == 0)
			return null;

		if (assignment.length != (1 + n * n))    // assumes assignment is'nt null
			throw new IllegalArgumentException("invalid input solution length");

		return createSolution(assignment, map, n);
	}


	public static int mapIndex(int[][] map, int n, int index) {
		// returns the relevant var with true value

		int step = index / n;    // returns map's row
		int city = index % n;    // returns index in map's row

		return map[step][city];
	}


	public static int[] createSolution(boolean[] assignment, int[][] map, int n) {
		// returns the solution tour of cities

		int[] solution = new int[n];

		for (int i = 1; i < assignment.length; i = i + 1) {    // first is irrelevant

			if (assignment[i]) {

				int mapIndex = mapIndex(map, n, i - 1);  // minus irrelevant one
				solution[(i - 1) / n] = (mapIndex - 1) % n;    // max n * n divided by 'n', map starts from one

			}
		}

		return solution;
	}


	// Task 19
	public static int[] solve(boolean[][] flights) {

		if (flights == null | !isLegalInstance(flights))
			throw new IllegalArgumentException("flights array is illegal");

		int n = flights.length;
		int[][] map = createVarsMap(n);

		SATSolver.init(n*n);
		encode(flights, map);

		boolean[] assignment = SATSolver.getSolution();

		if (assignment == null)				// enough checking time
			throw new IllegalArgumentException("TIMEOUT");

		return decode(assignment, map);    // assumes assignment != null
	}


	// Task 20
	public static boolean solve2(boolean[][] flights) {

		if (flights == null | !isLegalInstance(flights))
			throw new IllegalArgumentException("flights array is illegal");

		int [] solutionA = solve(flights);

		int counter = 1; int checks = 0;  	// decided to try two other options

		while (checks < 2 & counter < 2) {

			int [] solutionB = solve(flights);

			if (!sameSolutions(solutionA, solutionB))
				counter = counter + 1;      // another legal solutions

			checks = checks + 1;
		}

		return counter == 2;
	}


	public static boolean sameSolutions(int [] solA, int [] solB) {
	// checks if two solutions are the same

		if (solA == null | solB == null)        // insures nullity before lengths
			throw new IllegalArgumentException("TIMEOUT");

		boolean isSame = (solA.length == solB.length);

		boolean found;
		for (int a = 0 ; a < solA.length & isSame; a = a + 1) {

			int value = solA[a]; found = false;

			for (int b = 0 ; b < solA.length & !found ; b = b + 1){

				if (solB[b] == value)           // each value in 'a' appears in 'b'
					found = true;

				else isSame = false;            // not same if not found

			}

			if (found) isSame = true;

		}

		return isSame;
	}

}