import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintWriter;
import java.text.SimpleDateFormat;
import java.util.*;
import java.util.concurrent.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public class Solver {

    static int CO_METHOD_RANDOM = 0;
    static int CO_METHOD_ONE_POINT = 1;
    static int CO_METHOD_TWO_POINT = 2;

    static int MIG_STRATEGY_RANDOM = 0;
    static int MIG_STRATEGY_BEST_N = 1;

    WindFarmLayoutEvaluator wfle;
    int num_scenario;
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
    int migration_strategy;

    int[] crossover_methods;
    int[] tournament_sizes;
    double[] mutation_probability;

    Date started_at;
    Date now;
    long since_start;

    String csv_filename;

    public Solver(WindFarmLayoutEvaluator evaluator, int num_scenario) {
        this(evaluator, num_scenario, 5, 10000, MIG_STRATEGY_RANDOM, 2);
    }

    public Solver(WindFarmLayoutEvaluator evaluator, int num_scenario, int num_generations_to_migrate, int num_generations, int migration_strategy, int num_to_migrate) {
        started_at = new Date();
        wfle = evaluator;
        this.num_scenario = num_scenario;
        rand = new Random();
        grid = new ArrayList<>();

        // set up any parameter here, e.g pop size, cross_rate etc.
        num_islands = 3;
        num_individuals_per_island = 20;  // change this to anything you want

        this.num_generations_to_migrate = num_generations_to_migrate;
        this.num_generations = num_generations;
        min_generation_for_migration = 10;
        this.migration_strategy = migration_strategy;
        this.num_individuals_to_migrate = num_to_migrate;

        crossover_methods = new int[]{ CO_METHOD_RANDOM, CO_METHOD_ONE_POINT, CO_METHOD_TWO_POINT };
        tournament_sizes = new int[]{ 2, 2, 2 };
        mutation_probability = new double[]{ 0.05, 0.05, 0.05 };

        csv_filename = new SimpleDateFormat("yyyyMMdd-HHmmss").format(started_at);
        List<String> parts = Stream.of("output", csv_filename).collect(Collectors.toList());
        csv_filename = String.join("-", parts) + ".csv";
        output("Will write to " + csv_filename);
    }

    private void output(String message) {
        now = new Date();
        since_start = now.getTime() - started_at.getTime();
        System.out.println("[" + since_start + "] " + message);
    }
    
    private void write_csv_line(String[] line) {
        String line_to_write = String.join(",", line);
        try {
            PrintWriter pw = new PrintWriter(new FileOutputStream(
                    new File(csv_filename),
                    true
            ));

            pw.println(line_to_write);
            pw.close();
        } catch (FileNotFoundException ex) {
            output(" --- ! --- Cannot open CSV file for output --- ! ---");
            output(line_to_write);
        }
    }
   
    public void run_cw() {
        started_at = new Date();
        output("Starting");

        String[] header_row = new String[num_islands + 2];
        header_row[0] = "generation";
        header_row[1] = "global_best";

        for (int i=0; i < num_islands; i++)
            header_row[i+2] = "island_" + i + "_best";

        write_csv_line(header_row);

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
        double initial_fitness = evaluate();

        output("Initial population evaluated");

        String[] initial_output = new String[num_islands + 2];
        initial_output[0] = "0";
        initial_output[1] = Double.toString(initial_fitness);

        for (int island=0; island<num_islands; island++)
            initial_output[island+2] = Double.toString(get_best_from_island(island));

        write_csv_line(initial_output);

        /**** PUT YOUR OPTIMISER CODE HERE ***********/

        for (int generation=0; generation<num_generations; generation++) {
            // add some code to evolve a solution
            output("Generation " + (generation + 1) + " starting");

            for (int island_n=0; island_n<num_islands; island_n++) {
                output("\tIsland: " + island_n);
                // Choose some parents
                boolean[] left_parent = individuals[island_n][select_parent(island_n)];
                output("\tGot Left Parent");
                boolean[] right_parent = individuals[island_n][select_parent(island_n)];
                output("\tGot Right Parent");

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

                output("\tChild Created");

                // Mutation
                for (int i=0; i<grid.size(); i++) {
                    if (rand.nextDouble() < mutation_probability[island_n])
                        child[i] = !child[i];
                }

                output("\tMutation Complete");

                // Replacement
                double worst_fitness = 0;
                int worst_index = -1;
                for (int i = 0; i<num_individuals_per_island; i++) {
                    if (fits[island_n][i] > worst_fitness) {
                        worst_fitness = fits[island_n][i];
                        worst_index = i;
                    }
                }
                output("\tReplacement Target: " + worst_index);

                double child_fitness = evaluate_individual(child);
                // Only insert new child into population if it is better than the worst of the population
                if (child_fitness < worst_fitness && worst_index != -1) {
                    fits[island_n][worst_index] = child_fitness;
                    individuals[island_n][worst_index] = child;
                }

                output("\tReplaced");
            }

            // Migration
            if (generation >= min_generation_for_migration && generation % num_generations_to_migrate == 0) {
                output("Generation " + (generation + 1) + " - Migrate");
                for (int island = 0; island < num_islands; island++) {
                    int target_island = (island == num_islands - 1) ? 0 : island + 1;
                    if (migration_strategy == MIG_STRATEGY_RANDOM) {
                        int the_chosen_one = rand.nextInt(num_individuals_per_island);

                        int destination = get_worst_index_from_island(target_island);
                        individuals[target_island][destination] = individuals[island][the_chosen_one];
                    } else if (migration_strategy == MIG_STRATEGY_BEST_N) {

                    }
                }
            }

            output("Generation " + (generation + 1) + " evaluation");
            // Evaluate at the end of the generation
            double best_global_fitness = evaluate();

            String[] generation_output = new String[num_islands + 2];
            generation_output[0] = Integer.toString(generation + 1);
            generation_output[1] = Double.toString(best_global_fitness);

            for (int island=0; island<num_islands; island++)
                generation_output[island+2] = Double.toString(get_best_from_island(island));

            write_csv_line(generation_output);
            output("Generation " + (generation + 1) + " complete");
        }
    }

    private int select_parent(int island_n) {
        int best_parent_index = -1;
        double best_parent_fitness = Double.MAX_VALUE;

        // run a basic tournament
        for (int i=0; i<tournament_sizes[island_n]; i++) {
            int this_parent_index = rand.nextInt(num_individuals_per_island);
            double this_parent_fitness = fits[island_n][this_parent_index];

            if (this_parent_fitness < best_parent_fitness) {
                best_parent_index = this_parent_index;
                best_parent_fitness = this_parent_fitness;
            }
        }


        return best_parent_index;
    }

    private int get_worst_index_from_island(int island) {
        int worst_index_in_island = -1;
        double worst_in_island = 0;

        for (int i=0; i<num_individuals_per_island; i++) {
            if (fits[island][i] > worst_in_island) {
                worst_index_in_island = i;
                worst_in_island = fits[island][i];
            }
        }

        return worst_index_in_island;
    }

    private double get_best_from_island(int island) {
        double best_in_island = Double.MAX_VALUE;
        for (int i=0; i<num_individuals_per_island; i++)
            if (fits[island][i] < best_in_island)
                best_in_island = fits[island][i];

        return best_in_island;
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
    private double evaluate() {
        double minfit = Double.MAX_VALUE;

        ExecutorService executor = Executors.newFixedThreadPool(16);
        Future<Double>[][] tasks = new Future[num_islands][num_individuals_per_island];

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

                final double[][] layout_final = layout;
                final int nSc = this.num_scenario;

                Callable<Double> task = () -> {
                    String threadName = Thread.currentThread().getName();
//                    output("Hello from thread " + threadName);

                    // As the provided evaluator is not thread safe
                    // we create an instance of it local to this thread
                    KusiakLayoutEvaluator localWfle = new KusiakLayoutEvaluator();
                    WindScenario sc = new WindScenario("./Scenarios/practice_"+nSc+".xml");
                    localWfle.initialize(sc);

                    double coe;
                    if (localWfle.checkConstraint(layout_final)) {
                        localWfle.evaluate(layout_final);
                        coe = localWfle.getEnergyCost();
                    } else {
                        coe = Double.MAX_VALUE;
                    }
                    return coe;
                };

//                output("Saving " + island + "-" + p);
                tasks[island][p] = executor.submit(task);

//                fits[island][p] = coe;
//                if (fits[island][p] < minfit) {
//                    minfit = fits[island][p];
//                }

            }
        }

        try {
            // to check we're finished, get the final value first
            // This call blocks until it is done.
//            output("Going to get the final one (" + (num_islands-1) + "-" + (num_individuals_per_island-1) + ")");
            tasks[num_islands - 1][num_individuals_per_island - 1].get();
//            output("Got the final one");

            // Now that's done, we can set our fits like we had been doing before
            for (int island=0; island<num_islands; island++) {
                for (int individual=0; individual<num_individuals_per_island; individual++) {
//                    output("Fitting " + island + "-" + individual);
                    Double result = tasks[island][individual].get();

                    fits[island][individual] = result;

                    if (result < minfit) {
                        minfit = result;
                    }
                }
            }
        } catch (InterruptedException e) {
            e.printStackTrace();
        } catch (ExecutionException e) {
            e.printStackTrace();
        }

        output("\tBest Fitness: " + minfit);
        return minfit;
    }
}
