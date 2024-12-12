// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "ModelFactory.hpp"
#include "FixedBoundsConstraintsModel.hpp"
#include "HomogeneousEqualityConstrainedModel.hpp"
#include "BoundRelaxedModel.hpp"
#include "options/Options.hpp"

namespace uno {
   // note: ownership of the pointer is transferred
   std::unique_ptr<Model> ModelFactory::reformulate(std::unique_ptr<Model> model, const Options& options) {
      if (options.get_string("inequality_handling_method") == "primal_dual_interior_point") {
         // move the fixed variables to the set of general constraints
         if (not model->get_fixed_variables().empty()) {
            model = std::make_unique<FixedBoundsConstraintsModel>(std::move(model), options);
         }
         // if an equality-constrained problem is required (e.g. interior points or AL), reformulate the model with slacks
         model = std::make_unique<HomogeneousEqualityConstrainedModel>(std::move(model));
         // slightly relax the bound constraints
         model = std::make_unique<BoundRelaxedModel>(std::move(model), options);
      }
      return model;
   }
} // namespace