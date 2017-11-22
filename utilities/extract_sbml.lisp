;;/home/ubuntu/RiboCensus//libs/python_scripts/RiboCensus_annotate_otus.py
;;  % pathway-tools -lisp
;; (load "~/examples")
;; list the pgdb names in "list.txt"
(defun get-sbml(org)

(select-organism :org-id org)
(setq *sbml-export-file-name* (concatenate 'string org ".sbml")  )

(let* 
         (
                  (maybe-with-enzymes 
                       ( when (eq *sbml-reaction-subset* :only-with-enzymes)
                            (all-rxns :enzyme)
                       )
                   )

                  (initial-reaction-set (all-rxns *sbml-reaction-type*))
 
                  (reaction-set   initial-reaction-set)
         )

         ;; ... and export that set via the call to "2cobra"
          (  handler-case

            (2cobra *sbml-export-file-name* reaction-set)

             (file-error (c)   
                (clim-warn-user (format nil "Could not save to file ~A.~%Error: ~A"  *sbml-export-file-name* c))
             )
         )
)

)


(defun get-file (filename)
  (with-open-file (stream filename)
    (loop for line = (read-line stream nil)
          while line
          collect line)))

    ;;(message "sample ~A" n)

(dolist (n (get-file "list.txt")    )  
    (get-sbml  n)
)
