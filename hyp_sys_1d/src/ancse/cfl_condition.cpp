#include <ancse/cfl_condition.hpp>

#include <Eigen/Dense>


StandardCFLCondition <FVM>
:: StandardCFLCondition(const Grid &grid,
                        const std::shared_ptr<Model> &model,
                        double cfl_number)
    : grid(grid), model(model), cfl_number(cfl_number) {}

double StandardCFLCondition <FVM>
:: operator()(const Eigen::MatrixXd &u) const {

    auto n_cells = grid.n_cells;
    auto n_ghost = grid.n_ghost;


    Eigen::VectorXd lamdamax(n_cells-2*n_ghost);

    for (int i = n_ghost; i < n_cells-n_ghost;i++)
    {

        lamdamax(i-n_ghost) = model->max_eigenvalue(u.col(i));
    }

    double a_max = 0.0;

    for (int i = 0; i < n_cells-2*n_ghost; i++)
    {
        if (lamdamax(i)>a_max){
            a_max = lamdamax(i);
        }
    }
    //double a_max = lamdamax.maxCoeff();
    auto dt=cfl_number*grid.dx/a_max;

    return dt;

}

StandardCFLCondition <DG>
:: StandardCFLCondition(const Grid &grid,
                        const std::shared_ptr<Model> &model,
                        const DGHandler &dg_handler,
                        double cfl_number)
    : grid(grid), model(model), dg_handler (dg_handler),
      cfl_number(cfl_number) {}

double StandardCFLCondition <DG>
:: operator()(const Eigen::MatrixXd &u) const {

    auto n_cells = grid.n_cells;
    auto n_ghost = grid.n_ghost;

    return 0.0;
}


/// make CFL condition for FVM
std::shared_ptr<CFLCondition>
make_cfl_condition(const Grid &grid,
                   const std::shared_ptr<Model> &model,
                   double cfl_number) {

    return std::make_shared<StandardCFLCondition <FVM>>(grid, model,
                                                        cfl_number);

}

/// make CFL condition for DG
std::shared_ptr<CFLCondition>
make_cfl_condition(const Grid &grid,
                   const std::shared_ptr<Model> &model,
                   const DGHandler &dg_handler,
                   double cfl_number) {

    return std::make_shared<StandardCFLCondition <DG>>
            (grid, model, dg_handler, cfl_number);
}
