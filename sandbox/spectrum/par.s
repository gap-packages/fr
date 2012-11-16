	.file	"par.c"
	.section	.rodata.str1.1,"aMS",@progbits,1
.LC0:
	.string	"Hello World from thread %d\n"
.LC1:
	.string	"There are %d threads\n"
	.text
	.p2align 4,,15
	.type	main.omp_fn.0, @function
main.omp_fn.0:
	pushl	%ebp
	movl	%esp, %ebp
	pushl	%ebx
	subl	$20, %esp
	call	omp_get_thread_num
	movl	$.LC0, 4(%esp)
	movl	$1, (%esp)
	movl	%eax, %ebx
	movl	%eax, 8(%esp)
	call	__printf_chk
	call	GOMP_barrier
	testl	%ebx, %ebx
	je	.L6
	addl	$20, %esp
	popl	%ebx
	popl	%ebp
	.p2align 4,,2
	ret
	.p2align 4,,7
	.p2align 3
.L6:
	.p2align 4,,6
	call	omp_get_num_threads
	movl	8(%ebp), %edx
	movl	%eax, (%edx)
	movl	%eax, 8(%esp)
	movl	$.LC1, 4(%esp)
	movl	$1, (%esp)
	call	__printf_chk
	addl	$20, %esp
	popl	%ebx
	popl	%ebp
	ret
	.size	main.omp_fn.0, .-main.omp_fn.0
	.p2align 4,,15
.globl main
	.type	main, @function
main:
	pushl	%ebp
	movl	%esp, %ebp
	andl	$-16, %esp
	pushl	%ebx
	subl	$44, %esp
	leal	28(%esp), %ebx
	movl	%ebx, 4(%esp)
	movl	$0, 28(%esp)
	movl	$0, 8(%esp)
	movl	$main.omp_fn.0, (%esp)
	call	GOMP_parallel_start
	movl	%ebx, (%esp)
	call	main.omp_fn.0
	call	GOMP_parallel_end
	addl	$44, %esp
	xorl	%eax, %eax
	popl	%ebx
	movl	%ebp, %esp
	popl	%ebp
	ret
	.size	main, .-main
	.ident	"GCC: (Ubuntu 4.4.3-4ubuntu5) 4.4.3"
	.section	.note.GNU-stack,"",@progbits
