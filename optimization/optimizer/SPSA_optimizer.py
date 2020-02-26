from math import *
from pprint import pprint
from typing import (
    List,
    Optional,
    Tuple,
    )

import numpy as np

# sys.path.append('../examples/flatirons')
# import func_tools
# matplotlib.use('tkagg')
from optimization.data_logging.data_recorder import DataRecorder
from optimization.optimizer.ask_tell_optimizer import AskTellOptimizer
from optimization.optimizer.dimension.dimension_info import DimensionInfo


# import shapely

class SPSADimensionInfo:
    
    def __init__(self, theta: float, a_scale: float, info: DimensionInfo):
        self.theta: float = theta
        self.a_scale: float = a_scale
        self.info: DimensionInfo = info


class SPSAOptimizer(AskTellOptimizer):
    """
    A prototype implementation of a simultaneous perturbation stochastic approximation optimizer
    see https://www.jhuapl.edu/SPSA/
    see https://www.jhuapl.edu/SPSA/PDF-SPSA/Spall_An_Overview.PDF
    see https://en.wikipedia.org/wiki/Simultaneous_perturbation_stochastic_approximation
    see https://www.jhuapl.edu/SPSA/PDF-SPSA/Spall_Implementation_of_the_Simultaneous.PDF
    """
    
    def __init__(self,
                 a: float,  # s.t. a / (A+1)^alpha ~= smallest desired change in elements of theta in early iterations
                 c: float = 1e-3,  # ~= std dev of measurement noise
                 A: float = 1e2,  # <= # iterations expected, typically 10%
                 alpha: float = .602,  # or 1.0
                 gamma: float = .101,  # or 1/6th
                 dimensions: Optional[List[DimensionInfo]] = None,
                 num_estimates: int = 1,
                 ) -> None:
        self._recorder: Optional[DataRecorder] = None
        
        self._a = a
        self._c = c
        self._A = A
        self._alpha = alpha
        self._gamma = gamma
        
        self._k: int = 1  # step counter
        self._theta: Optional[np.ndarray] = None
        self._a_scale: Optional[np.ndarray] = None
        self._dimensions: Optional[List[DimensionInfo]] = dimensions
        self._num_estimates: int = num_estimates
        
        self._ck: float = self._c
    
    def setup(self, dimensions: [SPSADimensionInfo], recorder: DataRecorder) -> None:
        self._k = 1
        self._theta = np.fromiter((t.theta for t in dimensions), float)
        self._a_scale = np.fromiter((t.a_scale for t in dimensions), float)
        self._dimensions = [t.info for t in dimensions]
        self._ck: float = self._c
        
        self._recorder = recorder
        self._recorder.add_columns('generation', 'gradient_estimate')
    
    def ask(self, num: Optional[int] = None) -> [any]:
        num_estimates = ceil(num / 2) if num is not None else self._num_estimates
        num_dimensions = self.get_num_dimensions()
        ck = self._c / (self._k ** self._gamma)
        self._ck = ck
        result = []
        for c in range(num_estimates):
            positive_candidate = np.empty(num_dimensions)
            negative_candidate = np.empty(num_dimensions)
            for i in range(num_dimensions):
                delta_i = self._dimensions[i].sample()
                perturbation = ck * delta_i
                positive_candidate[i] = self._theta[i] + perturbation
                negative_candidate[i] = self._theta[i] - perturbation
            
            result.append(positive_candidate)
            result.append(negative_candidate)
        pprint(result)
        return result
    
    def tell(self, evaluations: [Tuple[float, any]]) -> None:
        gradient_estimate = np.zeros(self.get_num_dimensions())
        
        def convert_candidate(candidate):
            if isinstance(candidate, np.ndarray):
                return candidate
            return np.fromiter(candidate, dtype='float64')
        
        num_estimates = floor(len(evaluations) / 2)
        for i in range(num_estimates):
            base = i * 2
            positive_candidate = evaluations[base]
            negative_candidate = evaluations[base + 1]
            
            # compute gradient estimate
            
            
            performance_difference = positive_candidate[0] - negative_candidate[0]
            difference = convert_candidate(positive_candidate[1]) - convert_candidate(negative_candidate[1])  # = 2*ck*delta
            delta = performance_difference / (num_estimates * difference)
            if fabs(performance_difference) > 1e-12:
                gradient_estimate += delta
        
        ak = self._a / (self._k + self._A) ** self._alpha
        # print(ak * self.a_scale)
        update = ak * self._a_scale * gradient_estimate
        self._theta += update
        self._k += 1
        
        self._recorder.accumulate(evaluations, gradient_estimate)
    
    def best_solution(self) -> (Optional[float], any):
        return None, self._theta
    
    def get_num_candidates(self) -> int:
        return self._num_estimates * 2
    
    def get_candidate_block_size(self) -> int:
        return 2
    
    def get_num_dimensions(self) -> int:
        return len(self._dimensions)
