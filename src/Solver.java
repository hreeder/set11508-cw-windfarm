import java.util.ArrayList;
import java.util.Date;
import java.util.List;
import java.util.Random;

public class Solver {

    static int CO_METHOD_RANDOM = 0;
    static int CO_METHOD_ONE_POINT = 1;
    static int CO_METHOD_TWO_POINT = 2;

    WindFarmLayoutEvaluator wfle;
    boolean[][][] individuals;
    double[][] fits;
    Random rand;
    int num_islands;
    int num_individuals_per_island;
    
    ArrayList<double[]> grid;

    int num_generations_to_migrate;
    int num_individuals_to_migrate;
    int num_generations;
    int min_generation_for_migration;

    int[] crossover_methods;
    int[] tournament_sizes;

    Date started_at = new Date();
    Date now;
    long since_start;

    public Solver(WindFarmLayoutEvaluator evaluator) {
        wfle = evaluator;
        rand = new Random();
        grid = new ArrayList<double[]>();
        
        // set up any parameter here, e.g pop size, cross_rate etc.
        num_islands = 3;
        num_individuals_per_island = 10;  // change this to anything you want

        num_generations_to_migrate = 5;
        num_individuals_to_migrate = 2;
        num_generations = 10000;
        min_generation_for_migration = 10;

        crossover_methods = new int[num_islands];
        crossover_methods[0] = CO_METHOD_RANDOM;
        crossover_methods[1] = CO_METHOD_ONE_POINT;
        crossover_methods[2] = CO_METHOD_TWO_POINT;

        tournament_sizes = new int[num_islands];
        tournament_sizes[0] = 2;
        tournament_sizes[1] = 2;
        tournament_sizes[2] = 2;
    }

    private void output(String message) {
        now = new Date();
        since_start = now.getTime() - started_at.getTime();
        System.out.println("[" + since_start + "] " + message);
    }
    
    
   
    public void run_cw() {
        started_at = new Date();
        output("Starting");

        /************set up grid for scenario chosen  ***************************/
        // do not change or remove this section of code
        // centers must be > 8*R apart and avoid obstacles

        double interval = 8.001 * wfle.getTurbineRadius();

        for (double x=0.0; x<wfle.getFarmWidth(); x+=interval) {
            for (double y=0.0; y<wfle.getFarmHeight(); y+=interval) {
                boolean valid = true;
                for (int o=0; o<wfle.getObstacles().length; o++) {
                    double[] obs = wfle.getObstacles()[o];
                    if (x>obs[0] && y>obs[1] && x<obs[2] && y<obs[3]) {
                        valid = false;
                    }
                }

                if (valid) {
                    double[] point = {x, y};
                    grid.add(point);
                }
            }
        }

        /************initialize a population:*****************************/

        //  the variable grid.size() denotes the
        // maximum number of turbines for the given scenario

        individuals = new boolean[num_islands][num_individuals_per_island][grid.size()];
        fits = new double[num_islands][num_individuals_per_island];

        for (int island=0; island < num_islands; island ++) {
            for (int p=0; p < num_individuals_per_island; p++) {
                for (int i=0; i < grid.size(); i++) {
                    individuals[island][p][i] = rand.nextBoolean();
                }
            }
        }

        output("Population generated");

       /****** evaluate initial population  *************************/

        // this populates the fit[] array with the fitness values for each solution
        evaluate();

        output("Initial population evaluated");

        /**** PUT YOUR OPTIMISER CODE HERE ***********/

        for (int generation=0; generation<num_generations; generation++) {
            // add some code to evolve a solution

            for (int island_n=0; island_n<num_islands; island_n++) {
                // Choose some parents
                boolean[] left_parent = individuals[island_n][select_parent(island_n)];
                boolean[] right_parent = individuals[island_n][select_parent(island_n)];

                // Crossover
                boolean[] child = new boolean[grid.size()];

                if (crossover_methods[island_n] == CO_METHOD_RANDOM) {
                    for (int i=0; i < grid.size(); i++) {
                        boolean useRight = rand.nextBoolean();
                        child[i] = useRight ? right_parent[i] : left_parent[i];
                    }
                } else if (crossover_methods[island_n] == CO_METHOD_ONE_POINT) {
                    int midpoint = rand.nextInt(grid.size());
                    for (int i=0; i < grid.size(); i++) {
                        if (i < midpoint) {
                            child[i] = left_parent[i];
                        } else {
                            child[i] = right_parent[i];
                        }
                    }
                } else if (crossover_methods[island_n] == CO_METHOD_TWO_POINT) {
                    int left_midpoint = rand.nextInt(grid.size() / 2);
                    int right_midpoint = (grid.size() / 2) + rand.nextInt(grid.size() / 2);

                    for (int i=0; i < grid.size(); i++) {
                        if (i < left_midpoint || i >= right_midpoint) {
                            child[i] = left_parent[i];
                        } else {
                            child[i] = right_parent[i];
                        }
                    }
                }

                // Mutation
                // No mutation just now

                // Replacement
                double worst_fitness = 0;
                int worst_index = -1;
                for (int i = 0; i<num_individuals_per_island; i++) {
                    if (fits[island_n][i] > worst_fitness) {
                        worst_fitness = fits[island_n][i];
                        worst_index = i;
                    }
                }

                double child_fitness = evaluate_individual(child);
                // Only insert new child into population if it is better than the worst of the population
                if (child_fitness < worst_fitness && worst_index != -1) {
                    individuals[island_n][worst_index] = child;
                }
            }

            if (generation >= min_generation_for_migration && generation % num_generations_to_migrate == 0) {
                output("Generation " + generation + " - Migrate");
            }

            // Evaluate at the end of the generation
            evaluate();
            output("Generation " + generation + " complete");
        }
    }

    private int select_parent(int island_n) {
        int best_parent_index = -1;
        double best_parent_fitness = Double.MAX_VALUE;

        // run a basic tournament
        for (int i=0; i<tournament_sizes[island_n]; i++) {
            int this_parent_index = rand.nextInt(num_individuals_per_island);
            double this_parent_fitness = evaluate_individual(individuals[island_n][this_parent_index]);

            if (this_parent_fitness < best_parent_fitness) {
                best_parent_index = this_parent_index;
                best_parent_fitness = this_parent_fitness;
            }
        }


        return best_parent_index;
    }
    
    // evaluate a single chromosome
    private double evaluate_individual(boolean[] child) {
         int nturbines=0;
         for (int i=0; i<grid.size(); i++) {
                if (child[i]) {
                    nturbines++;
                }
         }

        double[][] layout = new double[nturbines][2];
        int l_i = 0;
        for (int i=0; i<grid.size(); i++) {
            if (child[i]) {
                layout[l_i][0] = grid.get(i)[0];
                layout[l_i][1] = grid.get(i)[1];
                l_i++;
            }
        }
	    
	    double coe;
	    if (wfle.checkConstraint(layout)) {
            wfle.evaluate(layout);
            coe = wfle.getEnergyCost();
            //System.out.println("layout valid");
	    } else {
		    coe = Double.MAX_VALUE;
	    }

     return coe;
	  
        
    }

    // evaluates the whole population
    private void evaluate() {
        double minfit = Double.MAX_VALUE;

        for (int island=0; island<num_islands; island++) {
            for (int p = 0; p < num_individuals_per_island; p++) {
                int nturbines = 0;
                for (int i = 0; i < grid.size(); i++) {
                    if (individuals[island][p][i]) {
                        nturbines++;
                    }
                }
    
                double[][] layout = new double[nturbines][2];
                int l_i = 0;
                for (int i = 0; i < grid.size(); i++) {
                    if (individuals[island][p][i]) {
                        layout[l_i][0] = grid.get(i)[0];
                        layout[l_i][1] = grid.get(i)[1];
                        l_i++;
                    }
                }
    
                double coe;
                if (wfle.checkConstraint(layout)) {
                    wfle.evaluate(layout);
                    coe = wfle.getEnergyCost();
                } else {
                    coe = Double.MAX_VALUE;
                }
    
                fits[island][p] = coe;
                if (fits[island][p] < minfit) {
                    minfit = fits[island][p];
                }
    
            }
        }
        output("\tBest Fitness: " + minfit);
    }

    
    
}
