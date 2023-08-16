# RocketOpt
Design optimization of a simplified rocket thrust chamber. This was a final project for Aero 588 Multidisciplinary Design Optimization taught by Professor Martins at the University of Michigan.

Gradient-based design optimization is applied to the design of a rocket thrust chamber with incrementally more constraints added. See DuvallJames_FinalProject.pdf for the final report which describes all of the sub-models and constraint scenarios and results in detail. 

See run_optimizer.ipynb for a condensed presentation and implementations for Isp-only and thrust-only constraint scenarios. Examples with more constraints following the final report to be added.

<!-- Contents
----------------

The repo should contain the following files:  

-----------------------------------
    RocketOpt
    ├── dataloader_poisson.py
    ├── dense_networks.py
    ├── fourier_layers.py    
    ├── hypernet_oneshot_common.py     
    ├── hypernet_oneshot_networks.py
    ├── hypernet_oneshot_train_poisson_dataset.py  
    ├── learning_rate_schedules.py
    ├── problem_settings.py
    ├── tf_util.py
    ├── train_hypernet.py
    ├── train_util.py
    └── README.md
----------------------------------- -->