#!/usr/bin/env python3
from fastq_read_stats import *
from multiprocessing import Process, Queue


def parse_seq_and_get_quality_stats(input_filepath: str) -> None:
    """
    Takes in a path to a fastq.gz file
    Generates quality stats and writes them directly
    to the stdout in tab-delimited format.
    """
    recs = pyfastx.Fastq(input_filepath, build_index=False)
    base_filename = os.path.basename(input_filepath)

    ## figure out the number of CPUs available to us
    NUMBER_OF_PROCESSING_PROCESSES = len(os.sched_getaffinity(0))

    ## set up queues for multiprocessing
    recordQueue = Queue()
    writerQueue = Queue()

    def fastq_to_queue(targetQueue):
        """
        Adds pyfastx objects for each read
        to the multiprocessing queue.

        Writes "DONE" when it's done so downstream
        processes can finish.

        Need a DONE for every CPU.
        """
        for r in recs:
            targetQueue.put(r)

        ## really make sure it finishes
        for _ in range(NUMBER_OF_PROCESSING_PROCESSES * 10):
            targetQueue.put("DONE")

    def process_queue_record(sourceQueue, targetQueue):
        """
        Reads pyfastx objects from the source queue
        Processes them and pushes the output to
        the target queue os we can write them out
        to the stdout.

        Writes "DONE" when its' done so downstream
        processes can finish.
        """
        while True:
            msg = sourceQueue.get()
            if msg == "DONE":
                break
            else:
                targetQueue.put(get_quality_stats(msg))
        targetQueue.put("DONE")

    def write_out_record(sourceQueue):
        """
        Reads quality stats objects from the source
        queue and writes TSV rows the stdout.
        """
        isHeaderWritten = False

        ## write out records
        while True:
            msg = sourceQueue.get()
            if msg == "DONE":
                break
            else:
                ## write header
                if not isHeaderWritten:
                    sys.stdout.write("\t".join(["input_filename", *msg._fields]))
                    sys.stdout.write("\n")
                    isHeaderWritten = True
                ## write rows
                sys.stdout.write("\t".join([base_filename, *[str(y) for y in msg]]))
                sys.stdout.write("\n")

    ## start process for reading fastq into the queue
    process_fastq_to_queue = Process(target=fastq_to_queue, args=[recordQueue])
    process_fastq_to_queue.daemon = True
    process_fastq_to_queue.start()

    ## start the processes for processing
    ## fastq objects from the queue.
    ## one process per CPU
    processing_processes_array = []

    for _ in range(NUMBER_OF_PROCESSING_PROCESSES):
        processing_processes_array.append(
            Process(target=process_queue_record, args=[recordQueue, writerQueue])
        )
    for p in processing_processes_array:
        p.daemon = True
        p.start()

    ## start the process for writing to
    ## the stdout
    process_write_out_record = Process(target=write_out_record, args=[writerQueue])
    process_write_out_record.daemon = True
    process_write_out_record.start()

    ## collect all the processes once they are done
    process_fastq_to_queue.join()

    for p in processing_processes_array:
        p.join()

    process_write_out_record.join()


def init_argparser():
    parser = argparse.ArgumentParser(
        prog="fastq_read_statsMP",
        description="Generates basic statistics given a fastq.gz as input",
    )

    parser.add_argument("input", help="Path to an input fastq.gz file")

    return parser


def main():
    args = init_argparser().parse_args()

    parse_seq_and_get_quality_stats(args.input)


if __name__ == "__main__":
    main()
