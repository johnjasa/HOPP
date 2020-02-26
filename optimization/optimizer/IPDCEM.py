import random
from typing import (
    Optional,
    Tuple,
    )

import numpy as np

# import shapely
from optimization.optimizer.DCEM_optimizer import DCEMOptimizer


# sys.path.append('../examples/flatirons')
# import func_tools
# matplotlib.use('tkagg')


class IPDCEM(DCEMOptimizer):
    """
    A prototype implementation of an incremental decomposed cross-entropy method.
    """
    
    def __init__(self, generation_size: int, selection_size: int, scale: float,
                 **kwargs
                 ) -> None:
        super().__init__(generation_size, 1.0, **kwargs)
        self._selection_size: int = selection_size
        self._scale: float = scale
        self._population: [Tuple[float, any]] = []
    
    def ask(self, num: Optional[int] = None) -> [any]:
        
        if len(self._population) == 0:
            return super().ask(num)
        
        if num is None:
            num = self._generation_size
        
        population = []
        for _ in range(num):
            base = self._population[random.randrange(len(self._population))][1]
            # candidate = [0.0] * len(self.dimensions)
            candidate = np.empty(self.get_num_dimensions())
            for i, dimension in enumerate(self._dimensions):
                candidate[i] = base[i] + (dimension.sample() - dimension.best_solution()) * self._scale
            population.append(candidate)
        
        return population
    
    def tell(self, evaluations: [Tuple[float, any]]) -> None:
        print('eval: ', [sample[0] for sample in evaluations])
        self._population.extend(evaluations)
        self._population.sort(key=lambda evaluation: evaluation[0], reverse=True)
        # self.population = self.population[0:self.selection_size]
        del self._population[self._selection_size:]
        print('pop: ', [sample[0] for sample in self._population])
        
        for i, dimension in enumerate(self._dimensions):
            dimension.update([evaluation[1][i] for evaluation in self._population])
    
    def best_solution(self) -> (Optional[float], any):
        return self._population[0] if len(self._population) > 0 else super().best_solution()
