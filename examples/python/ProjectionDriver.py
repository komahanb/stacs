from pspace import PSPACE

class ProjectionDriver:
    def __init__(self, parameter_container, callback_execute_solver, callback_compute_moments):
        self.parameter_container = parameter_container
        self.execute_solver      = callback_execute_solver
        self.compute_moments     = callback_compute_moments
    def execute(self):
        self.solver = self.execute_solver(None, self.parameter_container)
    def moments(self):
        return self.compute_moments(self.solver, self.solver.getNumTimeSteps() + 1, self.parameter_container.getNumBasisTerms())
    
class SamplingDriver:
    def __init__(self, parameter_container, callback_execute_solver):
        self.parameter_container = parameter_container
        self.execute_solver      = callback_execute_solver
    def execute(self):
        Q = self.parameter_container.getNumQuadraturePoints()
        for q in range(Q):
            wq, zq, yq = self.parameter_container.quadrature(q)
            self.solver = self.execute_solver(yq, None)
