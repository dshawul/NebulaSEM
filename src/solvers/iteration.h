#ifndef __ITERATION_H
#define __ITERATION_H

/**
  Iteration object that does common book keeping stuff
  for all solvers.
 */
class Iteration {
    private:
        Int starti; /**< start index */
        Int endi; /**< end index */
        Int i; /**< current iteration index */
        Int n_deferred; /**< number of deferred correction steps */
        Int idf; /**< current deferred correction index */
        bool skip; /**< skip iteration */

    public:
        Iteration(int step) {
            skip = false;
            if(step < -Constants::EqualEpsilon) {
                skip = true;
                step = 0;
            }
            starti = Controls::write_interval * step + 1;
            endi = Controls::write_interval * (step + Controls::amr_step);
            if(endi > Controls::end_step) endi = Controls::end_step;
            n_deferred = Controls::n_deferred;
            i = starti;
            Controls::current_step = i;
            idf = 0;
            if(MP::printOn)
                std::cout << "--------------------------------------------\n";
            Mesh::read_fields(step);
            Mesh::getProbeCells(Mesh::probeCells);
            forEachCellField (initTimeSeries());
            if(MP::printOn) {
                std::cout << "--------------------------------------------\n";
                MP::printH("Starting iterations.\n");
            }
        }
        bool start() {
            return (i == starti);
        }
        bool end() {
            /*zero iterations*/
            if(skip) {
                Mesh::write_fields(0);
                return true;
            }
            /*check*/
            if(i > endi)
                return true;
            /*iteration number*/
            if(MP::printOn && idf == 0) {
                if(Controls::state == Controls::STEADY)
                    MP::printH("Step %d\n",i);
                else
                    MP::printH("Time %f\n",i * Controls::dt);
            }
            return false;
        }
        void next() {
            idf++;
            if(idf <= n_deferred)
                return;
            idf = 0;

            /*set output printing*/
            MP::printOn = (MP::host_id == 0 && 
                    (MP::hasElapsed(Controls::print_time) || i == Controls::end_step - 1)); 

            /*update time series*/
            forEachCellField(updateTimeSeries(i));

            /*write result to file*/
            if((i % Controls::write_interval) == 0) {
                Int step = i / Controls::write_interval;
                Mesh::write_fields(step);
            }

            /*increment*/
            i++;
            Controls::current_step = i;
        }
        ~Iteration() {
        }
        Int get_step() {
            return i;
        }
};
/**
  Iterator for AMR
 */
class AmrIteration {
    private:
        int starti; /**< start index */
        int endi; /**< end index */
        int i; /**< current index */

    public:
        AmrIteration() {
            starti = Controls::start_step / Controls::write_interval - Controls::amr_step;
            endi = Controls::end_step / Controls::write_interval;
            if(endi < 1) endi = 1;
            if(!Controls::amr_step)
                Controls::amr_step = endi;
            i = starti;

            int in = (i < 0) ? 0 : i;
            if (MP::n_hosts > 1)
                Prepare::decomposeMesh(in);
            else
                Mesh::LoadMesh(in);
        }
        bool start() {
            return (i == starti);
        }
        bool end() {
            if(i >= endi)
                return true;
            return false;
        }
        void next() {
            i += Controls::amr_step;
            if(i < endi) {
                if (MP::host_id == 0) {
                    if(MP::n_hosts > 1) {
                        System::cd(MP::workingDir);
                        Prepare::mergeFields(i);
                    }

                    Prepare::refineMesh(i);
                }

                if(MP::n_hosts > 1)
                    Prepare::decomposeMesh(i);
                else
                    Mesh::LoadMesh(i);  
            }
        }
        ~AmrIteration() {
        }
        int get_step() {
            return i;
        }
};

#endif
