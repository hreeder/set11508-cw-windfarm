import Solver.*

data class SolverSettings(
        val evaluator: WindFarmLayoutEvaluator,
        val scenario_filename: String = "./Scenarios/practice_0.xml",
        val max_hours_to_run: Int = 8,
        val max_generations: Int = 10000,
        val num_generations_to_migrate: Int = 5,
        val migration_strategy: Int = MIG_STRATEGY_RANDOM,
        val num_to_migrate: Int = 0,
        val min_generation_for_migration: Int = 10,

        val crossover_methods: Array<Int> = arrayOf(CO_METHOD_RANDOM, CO_METHOD_ONE_POINT, CO_METHOD_TWO_POINT),
        val tournament_sizes: Array<Int> = arrayOf(2, 2, 2),
        val mutation_probability: Array<Double> = arrayOf(0.05, 0.05, 0.05)
)