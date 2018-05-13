# HPCfinalproject-psoriasis.c

1.Select an article from an academic journal that describes a Cellular Automaton model and submit your selection to me for approval by Tuesday, 24-April.

2.Write a serial version of the code that implements the model in either C or Python. You must certainly won’t be writing from scratch -- I strongly encourage you to take a previous example from the course and modify it as needed. 

3.Using either OpenMP (C) or MPI (C or Python) create a parallelized version of this code. No animations are necessary.  

4.The project is due no later than 5 pm on May 14th. There are no extensions. 

# This project is a simulation of how psoriasis generates. Depends on the article I choose, the development of psoriasis can be divided into three parts: beginning(small point), spreading(larger area), healing(new skin grow out from the middle, shape becomes circular ring). So, I wrote the code to simulate one psoriasis. Green area are good skin, red points is the bad skin spreading out the disease, grey skin is the skin being infected, and then new green skin came out at the middle when the psoriasis reached a size that the red points do not have the power to keep infecting other good skin. Eventually, the skin cured. But the skin uninfected and the skin new generated have some differeces. The different symbols differenciate them. So, that it! My final project.  
