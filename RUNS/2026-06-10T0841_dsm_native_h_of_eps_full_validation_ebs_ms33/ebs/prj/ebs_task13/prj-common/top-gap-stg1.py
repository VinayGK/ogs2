try:
    import ogs.callbacks as OpenGeoSys
except ModuleNotFoundError:
    import OpenGeoSys
#below code modified from Ogata Banks benchmark
class BCDirichletStg1(OpenGeoSys.BoundaryCondition):
    def getDirichletBCValue(self, t, coords, node_id, primary_vars):
        current_displacement = primary_vars[2]
        if current_displacement < 5.0e-3:
          #print('inside python BC1')
          return (False, 0.0)
        return (True, 5.0e-3)
# instantiate the BC objects used by OpenGeoSys
dirichlet_top1 = BCDirichletStg1()
